__name__         = "PanGIA-VIS: PanGIA result visualization tool: class-based Bokeh implementation"
__author__       = "August (Gus) Thomas"
__adaptedfrom__  = "PanGIA-VIS: PanGIA result visualization tool, By: Po-E (Paul) Li, Bioscience Division, Los Alamos National Laboratory"
__version__      = "0.8.2"
__date__         = "2021/06/15"

import pandas as pd
import numpy as np
import sys
import os
import json
from math import pi
from operator import itemgetter
import re
import hashlib
from os.path import dirname, join, isfile
from bokeh.io import curdoc
from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row, layout, widgetbox, column
from bokeh.models import ColumnDataSource, HoverTool, Div, FactorRange, Range1d, TapTool, CustomJS, ColorBar, CDSView, GroupFilter, MultiSelect
from bokeh.models.widgets import Button, Slider, Select, TextInput, RadioButtonGroup, DataTable, TableColumn, NumberFormatter, Panel, Tabs, HTMLTemplateFormatter
from bokeh.models.callbacks import CustomJS
from bokeh.palettes import RdYlBu11, Spectral4, Spectral11, Set1, Viridis256, Turbo256
from bokeh.util import logconfig
from bokeh.resources import INLINE
from bokeh.embed import components, autoload_static
from bokeh.transform import factor_cmap, linear_cmap
from app.models import Results, ResultsToItem, Item
from pathlib import Path
import math
from sklearn import preprocessing

# ----------------------------------------------------------------------------------------------------------------------
# BOKEH CLASS
# ----------------------------------------------------------------------------------------------------------------------
class BokehObject():

    ###############################################################################
    # init global variables --- not sure if we need all of these (probably not)
    ###############################################################################
    ranks = ['strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']

    tsv_md5sum     = ""            # result tsv md5sum
    refresh_period = 0   # refresh period in seconds
    banner_div     = Div(text = "", width = 800)  # default page

    COLORS         = RdYlBu11 #RdYlBu11
    GRADIENT       = Turbo256
    C_SET1         = Set1[9]
    MIN_SIZE       = 10
    MAX_SIZE       = 36

    GCOV_CACHE     = {}
    BSAT_CACHE     = {}
    MASK_CACHE     = {}

    bgfiles        = []
    bg_mask        = {}

    # An alternative to the config_setup() method. Worth investigating as another user-interaction,
    # governed by callbacks within the visualizer itself. Would also make this script slightly neater.
    config_file = join(dirname(__file__),"config.json")

    # Initialize individual object attributes:
    # NOTE: To make routes.py as neat as possible, the majority of methods in this class map to object attributes initialized here.
    def __init__(self, result):

        # Support & File-Loading Attributes: 
        self.result        = result
        self.filepath      = result.path
        self.configs       = self.config_setup()
        self.report_df     = self.fetch_report_df()
        self.log_text      = self.fetch_log_text()
        self.hover         = self.define_tooltips()

        # Data Source Attributes:
        self.source        = self.test_source()
        self.read_info     = self.parseLog() # Fetch relevant log info for piecharts

        # Figure-Specific Attributes:
        self.dotplot       = self.dotplot_skeleton()
        self.piechart_list = self.piechart_sources()
        self.datatable     = self.construct_dt()

    ###############################################################################
    # DataSource Construction - Used in all Figure Methods:
    ###############################################################################

    # This method creates self.source and fills out source.data, which are both used in many other places by this class.
    # Note that self.source and it's attribute source.data map directly to the .tsv file produced by PanGIA.
    def test_source(self):

        source                     =  ColumnDataSource()
        source.data                =  dict(

            
            # Name and Level:
            taxa                   = self.report_df["NAME"].tolist(),
            taxid                  = [ str(x).replace(".0","") for x in self.report_df["TAXID"].values ],
            rank_level             = self.report_df["LEVEL"].to_list(),

            # Counts:
            raw_rd_cnt             = self.report_df["READ_COUNT"].to_list(),
            tol_genome_sz          = self.report_df["TOL_GENOME_SIZE"].to_list(),
            relative_abun          = self.report_df["REL_ABUNDANCE"].to_list(), #ra: 'relative abundance'

            # Scoring:
            min_score              = self.report_df["SCORE"].to_list(),
            score_bg               = self.report_df["SCORE_BG"].tolist() if 'SCORE_BG' in self.report_df else [None]*len(c), #bg: 'background'
            score_uniq             = self.report_df["SCORE_UNIQ"].to_list(), #score_uniq: 'score based on phylogenic uniqueness'

            # Normalized Data:
            norm_rd_cnt            = self.report_df["READ_COUNT_RNR"].to_list(), #rnr: 'normalized read count'
            min_norm_rd_cnt_combo      = self.report_df["READ_COUNT_RSNB"].to_list(), #rsnb: 'normalized & combined (primary + secondary hits) read count'
            read_primary           = self.report_df["PRI_READ_COUNT"].to_list(), #pri: 'primary read'
            rpkm                   = self.report_df["RPKM"].to_list(), #rpkm: 'reads per kilobase per million mapped reads'

            # Coverage-Specific Data:
            linear_len             = self.report_df["LINEAR_LENGTH"].to_list(), # lnr_len: 'linear length'
            min_linear_cov         = self.report_df["LINEAR_COV"].to_list(), #lnr_cov: 'linear coverage'
            rankspec_min_depth_cov = self.report_df["RS_DEPTH_COV_NR"].to_list() if "RS_DEPTH_COV_NR" in self.report_df else [None]*len(c), #rsdcnr: 'normalized rank specific depth-of-coverage'                         
            min_depth_cov               = self.report_df["DEPTH_COV"].to_list(), #depth_cov: 'depth-of-coverage'

            # Rank-Specific Data:
            strain                 = self.report_df["STR"].to_list(),
            species                = self.report_df["SPE"].to_list(),
            genus                  = self.report_df["GEN"].to_list(),
            family                 = self.report_df["FAM"].to_list(),
            order                  = self.report_df["ORD"].to_list(),                                 
            clade                  = self.report_df["CLA"].to_list(),
            phylum                 = self.report_df["PHY"].to_list(),
            superkingdom           = self.report_df["SK"].to_list(),
            root                   = self.report_df["ROOT"].to_list(),


            # Normalized Rank-Specific Data:
            strain_rnr             = self.report_df['STR_rnr'].to_list(),
            species_rnr            = self.report_df['SPE_rnr'].to_list(),
            genus_rnr              = self.report_df['GEN_rnr'].to_list(),
            family_rnr             = self.report_df['FAM_rnr'].to_list(),
            order_rnr              = self.report_df['ORD_rnr'].to_list(),
            clade_rnr              = self.report_df['CLA_rnr'].to_list(),
            phylum_rnr             = self.report_df['PHY_rnr'].to_list(),
            superkingdom_rnr       = self.report_df['SK_rnr'].to_list(),
            root_rnr               = self.report_df['ROOT_rnr'].to_list(),

            # Normalized + Combined (Primary + Secondary) Rank-Specific Data:
            strain_rnb             = self.report_df['STR_rnb'].to_list(),
            species_rnb            = self.report_df['SPE_rnb'].to_list(),
            genus_rnb              = self.report_df['GEN_rnb'].to_list(),
            family_rnb             = self.report_df['FAM_rnb'].to_list(),
            order_rnb              = self.report_df['ORD_rnb'].to_list(),
            clade_rnb              = self.report_df['CLA_rnb'].to_list(),
            phylum_rnb             = self.report_df['PHY_rnb'].to_list(),
            superkingdom_rnb       = self.report_df['SK_rnb'].to_list(),
            root_rnb               = self.report_df['ROOT_rnb'].to_list(),

            # Rank-Specific Taxonomic Specificity Level Data:
            strain_ri              = self.report_df['STR_ri'].to_list()  if 'STR_ri'   in self.report_df else [None]*len(c),
            species_ri             = self.report_df['SPE_ri'].to_list()  if 'SPE_ri'   in self.report_df else [None]*len(c),
            genus_ri               = self.report_df['GEN_ri'].to_list()  if 'GEN_ri'   in self.report_df else [None]*len(c),
            family_ri              = self.report_df['FAM_ri'].to_list()  if 'FAM_ri'   in self.report_df else [None]*len(c),
            order_ri               = self.report_df['ORD_ri'].to_list()  if 'ORD_ri'   in self.report_df else [None]*len(c),
            clade_ri               = self.report_df['CLA_ri'].to_list()  if 'CLA_ri'   in self.report_df else [None]*len(c),
            phylum_ri              = self.report_df['PHY_ri'].to_list()  if 'PHY_ri'   in self.report_df else [None]*len(c),
            superkingdom_ri        = self.report_df['SK_ri'].to_list()   if 'SK_ri'    in self.report_df else [None]*len(c),
            root_ri                = self.report_df['ROOT_ri'].to_list() if 'ROOT_ri'  in self.report_df else [None]*len(c),

            # Pathogen-Specific Data:
            patho_host             = self.report_df["HOST"].to_list(), #p_host: name of a pathogen's host
#           patho_origin           = self.report_df["p_src"].to_list(), #p_src: source organism of given pathogen
            patho_location         = self.report_df["LOCATION"].to_list(), #p_loc: (genomic) location of given pathogen
            patho_disease          = self.report_df["DISEASE"].to_list(), #p_dse: disease(s) caused by given pathogen
            pathogen               = ["Yes" if x == "Pathogen" else "No" for x in self.report_df['PATHOGEN'].values]
            )

        norm_value = self.report_df["REL_ABUNDANCE"].to_list()
        minnorm, maxnorm = 0.1, 0.4

        for i, val in enumerate(norm_value):
            norm_value[i] = abs((val - minnorm)) / (maxnorm - minnorm)

        source.data["NORM_REL_ABUNDANCE"] = norm_value

        return source

    ###############################################################################
    # Figure Construction - Dotplot, Piecharts, and DataTable:
    ###############################################################################

    ######
    ###### CONSTRUCT DOTPLOT:
    ######

    # Constructs a dotplot centered in our webpage - the "main attraction" as it were.
    # Like our other figures, dotplot's behavior will be governed by a CustomJS callback upon user interaction or detection of new data.
    def dotplot_skeleton(self):

        source = self.source
        dot_hover = self.hover[0]

        # So color_map (and eventually size_map) need to be adjustable via a Button Control in the JS Callback.
        color_map = linear_cmap(field_name = 'min_score', palette = Turbo256, low = min(source.data['min_score']), high = max(source.data['min_score']))

#        size_limit=[10]*len(self.report_df['NAME'].values)
#        max_x = max(self.report_df["DEPTH_COV"].values) if max(self.report_df["DEPTH_COV"].values) > 0 else 0.001
#        size_limit = [int(self.MIN_SIZE + x/max_x*(self.MAX_SIZE - self.MIN_SIZE)) for x in self.report_df["DEPTH_COV"].values]

#        species_view = CDSView(source = self.source, filters = [GroupFilter(column_name = 'species', group='species')])

        dotplot = figure(plot_width = self.configs['displays']['dot_plot_width'],                     plot_height    = self.configs['displays']['dot_plot_height'], 
                        x_range = FactorRange(*self.report_df["NAME"].unique()), output_backend = self.configs['displays']['output_backend'],
                        tools       = ["pan, wheel_zoom, box_zoom, reset, tap", dot_hover]
                        )

        dotplot.circle(x = "taxa", y = "raw_rd_cnt", source = source, radius = "NORM_REL_ABUNDANCE", color = color_map, line_color = color_map,
                    fill_alpha = 0.6, line_alpha = 0.7)

        # To make color_bar adjustable, the title needs to be populated by the call to JSCallback.
        # See class methods js_link() and js_on_change() for help.
        
        color_bar = ColorBar(color_mapper = color_map['transform'], width = 8, label_standoff = 12, title = 'Score Gradient',
            title_text_align = 'left', title_text_font_size = '20px')
        dotplot.add_layout(color_bar, 'right')
        
        dotplot.xaxis.axis_label              = "Taxa"                             
        dotplot.xaxis.major_label_orientation = 3.1415926/4
        dotplot.yaxis.axis_label              = "Read Count"

        #toolbar_location="below", toolbar_sticky=False,

        return dotplot

    ######
    ###### CONSTRUCT PIECHARTS:
    ######

    def piechart_sources(self):

#       pieHover = self.tooltips[1]

        ######
        ###### CONSTRUCT PIECHART DATA SOURCES:
        ######

        ###### Pie Chart Data Source Build:
        # Because our piecharts are constructed differently from dotplot, we must specify both source and source.data here.
        # We opted to construct source.data for each by referencing the index of a list, called 'pievalue_list', returned by a method in the 'Supporting Methods' section. 
        pieInReadsDS = ColumnDataSource(data = dict( name = ['NA'], start_angle = [0], end_angle = [2*pi], color = ['#EFF0F1'], val = ['NA'], pct = ['NA']))
        pieFlagDS    = ColumnDataSource(data = dict( name = ['NA'], start_angle = [0], end_angle = [2*pi], color = ['#EFF0F1'], val = ['NA'], pct = ['NA']))
        piePathoDS   = ColumnDataSource(data = dict( name = ['NA'], start_angle = [0], end_angle = [2*pi], color = ['#EFF0F1'], val = ['NA'], pct = ['NA']))

        ###### In-Reads Data Source:

        logfile            = self.read_info
        pievalue_list      = self.genPieValues(logfile)
        pieInReadsDS.data  = dict(
                     name  = pievalue_list[0], start_angle = pievalue_list[1], end_angle = pievalue_list[2],
                     color = pievalue_list[3], val         = pievalue_list[4], pct       = pievalue_list[5])

        ###### Pathos Data Source:

        info                 = {}
        pievalue_list        = []
        info['Pathogen']     = self.report_df.loc[self.report_df.PATHOGEN == "Pathogen","PRI_READ_COUNT"].sum()
        info['Not pathogen'] = self.report_df.loc[:,"PRI_READ_COUNT"].sum() - info['Pathogen']
        pievalue_list        = self.genPieValues(info)
        piePathoDS.data      = dict(
            name  = pievalue_list[0], start_angle = pievalue_list[1], end_angle = pievalue_list[2],
            color = pievalue_list[3], val         = pievalue_list[4], pct       = pievalue_list[5])

        ###### Flag Data Source:

        if "FLAG" in self.report_df:
            info          = {}
            piechart_list = []
            for flag in self.report_df.FLAG.unique():
                if flag   == 'B':
                    name  = "Bacteria"
                elif flag == 'A':
                    name  = "Archae"
                elif flag == 'E':
                    name  = "Eukaryota"
                elif flag == 'V':
                    name  = "Viruses"
                else:
                    name  = flag
                
                info[name] = self.report_df.loc[self.report_df.FLAG == flag, "PRI_READ_COUNT"].sum()

            pievalue_list = self.genPieValues(info)
            
            pieFlagDS.data  = dict(
                      name  = pievalue_list[0], start_angle = pievalue_list[1], end_angle = pievalue_list[2],
                      color = pievalue_list[3], val         = pievalue_list[4], pct       = pievalue_list[5])

        ######
        ###### BUILD OUT PIECHARTS:
        ######

        ###### Log File In-Reads Pie Chart:
        pieInReadsFigure   = figure(
            x_range        = (-1.3, 4),
            output_backend = self.configs['displays']['output_backend'], 
            y_range        = (-2, 2), 
            plot_width     = self.configs['displays']['dashboard_pie_width'], 
            plot_height    = self.configs['displays']['dashboard_pie_height'], 
            title          = "Total reads:",
            tools          = [self.hover[1]]
            )

        pieInReadsFigure.annular_wedge(
            x = 0, y = 0, alpha = 0.7,
            legend_label = 'name', start_angle  = 'start_angle', end_angle = 'end_angle', color = 'color',
            inner_radius = 0.7,    outer_radius = 1.2,           source    = pieInReadsDS
            )

        ###### Flag Distribution Piechart:
        pieFlagFigure   = figure(
            x_range     = (-1.5, 4), 
            y_range     = (-2, 2), 
            plot_width  = self.configs['displays']['dashboard_pie_width'], 
            plot_height = self.configs['displays']['dashboard_pie_height'], 
            title       = "Target reads distribution:", 
            tools       = [self.hover[1]]
            )

        pieFlagFigure.annular_wedge(
            x = 0, y = 0, alpha = 0.7,
            legend_label = 'name', start_angle  = 'start_angle', end_angle = 'end_angle', color = 'color',
            inner_radius = 0.7,    outer_radius = 1.2,           source    = pieFlagDS
            )

        ###### Pathogen Stats Piechart:
        piePathoFigure      = figure( 
            x_range         = (-1.27, 4), 
            output_backend  = self.configs['displays']['output_backend'], 
            y_range         = (-2, 2), 
            plot_width      = self.configs['displays']['dashboard_pie_width'], 
            plot_height     = self.configs['displays']['dashboard_pie_height'], 
            title           = "Pathogen reads distribution:",
            tools           = [self.hover[1]]
            )

        piePathoFigure.annular_wedge(
            x = 0, y = 0, alpha = 0.7,
            legend_label = 'name', start_angle  = 'start_angle', end_angle = 'end_angle', color = 'color',
            inner_radius = 0.7,    outer_radius = 1.2,           source    = piePathoDS
            )

        ######
        ###### RETURN FORMATTED PIECHARTS:
        ######

        piechart_list = [pieInReadsFigure, pieFlagFigure, piePathoFigure]

        for chart in piechart_list:
            chart.axis.visible       = False
            chart.grid.visible       = False
            chart.legend.location    = "center_right"
            chart.toolbar.logo       = None
            chart.toolbar_location   = None
            chart.outline_line_width = 0
            chart.outline_line_alpha = 0

        return piechart_list

    ######
    ###### CONSTRUCT DATA TABLE:
    ######

    def construct_dt(self):
        table_cols = [
            TableColumn(field = "taxa",                   title = "Name", width = 800),
            TableColumn(field = "rank_level",             title = "Rank"),
            TableColumn(field = "norm_rd_cnt",            title = "Normalized Read Count",                       formatter = NumberFormatter(format = '0,0.00')),
            TableColumn(field = "min_norm_rd_cnt_combo",  title = "Normalized & Combined Read Count",            formatter = NumberFormatter(format = '0,0.00')),
            TableColumn(field = "read_percent_id",        title = "Read Percent Identity",                       formatter = NumberFormatter(format = '0,0')),
            TableColumn(field = "min_score",              title = "Score",                                       formatter = NumberFormatter(format = '0.00')),
            TableColumn(field = "score_uniq",             title = "Score (Unique)",                              formatter = NumberFormatter(format = '0.00')),
            TableColumn(field = "score_bg",               title = "Score (Background)",                          formatter = NumberFormatter(format = '0.00')),
            TableColumn(field = "rpkm",                   title = "(Reads/kb)/1M Mapped Reads",                  formatter = NumberFormatter(format = '0,0')),
            TableColumn(field = "min_linear_cov",         title = "Genome Coverage",                             formatter = NumberFormatter(format = '0,0.00')),
            TableColumn(field = "rankspec_min_depth_cov", title = "Normalized Rank-Specific Depth-of-Coverage",  formatter = NumberFormatter(format = '0,0.00')),
            TableColumn(field = "min_depth_cov",          title = "Depth-of-Coverage",                           formatter = NumberFormatter(format = '0,0.00')),
            TableColumn(field = "relative_abun",          title = "Relative Abundance",                          formatter = NumberFormatter(format = '0.00%')),
            TableColumn(field = "pathogen",               title = "Pathogen"),
            ]

        #use "index_position=None" instead of "row_headers=False" on bokeh ver >0.12.10
        data_table = DataTable(
                source         = self.source, 
                columns        = table_cols, 
                index_position = None, 
                width          = self.configs['displays']['total_width'],
                height         = 200
                )
        result_table = column(data_table)
        return result_table


    ###############################################################################
    # Custom JavaScript Callback Methods
    ###############################################################################

    # All JS Callbacks require a dictionary called "controls" - it isn't clear whether a single "controls" dictionary with everything
    # in it will do, or if we need a method that returns a different dictionary each time it is called, customized for each JS callback event.
    # In the example, the "new data" seemed ready to accept everything in source. 

    def control_change(self):
        
        # These variables are used for dropdown menus in the interface - maybe could also be used for box-checkers.
        axis_options = [
            'None',
            'Raw Read Count',
            'Raw Read Count: Primary Alignment Only',
            'Normalized Read Count (RNR)',
            'Normalized Rank-Specific Read Count (RSNB)',
            'Genome Coverage',
            'Reads per Kilobase per Million Mapped Reads (RPKM)',
            'Normalized Rank-Specific Depth-of-Coverage (RSNR)',
            'Depth-of-Coverage',
            'Score',
            'Unique Score',
            'Backgound Score',
            'Relative Abundance'
            ]

        y_axis_options = [
            'Raw Read Count',
            'Raw Read Count: Primary Alignment Only',
            'Normalized Read Count (RNR)',
            'Normalized Rank-Specific Read Count (RSNB)',
            'Reads per Kilobase per Million Mapped Reads (RPKM)',
            'Normalized Rank-Specific Depth-of-Coverage (RSNR)',
            'Depth-of-Coverage',
            'Relative Abundance'
            ]

        size_comparison_options = [
            'None',
            'Genome Coverage',
            'Score',
            'Normalized Rank-Specific Depth-of-Coverage (RSNR)',
            'Depth-of-Coverage',
            'Relative Abundance'
            ]

        color_conparison_options = [
            'None',
            'Genome Coverage',
            'Score',
            'Relative Abundance'
            ]

        rank_options = ["strain", "species", "genus", "family",
                        "order", "class", "phylum",  "superkingdom", "root"
                       ]

        # To this dict we should add ALL the things that can be changed in the viewer.
        controls = {
            "raw_rd_cnt": Slider(             title = "Minimum Raw Read Count",                   value = 0,  start = 0,   end = 500,   step = 1),
            "min_linear_cov": Slider(         title = "Minimum Linear Coverage",                  value = self.configs['cutoffs']['def_val_min_cov'],    start = 0,   end = 1,     step = 0.01),
            "min_score": Slider(              title = "Minimum Score",                            value = self.configs['cutoffs']['def_val_min_score'],  start = 0, end = 1,     step = 0.1),
            "min_norm_rd_cnt_combo": Slider(  title = "Minimum Normalized & Combined Read Count", value = self.configs['cutoffs']['def_val_max_r_rsnb'], start = 0,   end = 100,   step = 1),
            "min_depth_cov": Slider(          title = "Minimum Depth-of-Coverage",                value = self.configs['cutoffs']['def_val_min_dc'],     start = 0,   end = 10000, step = 10),
            "rankspec_min_depth_cov": Slider( title = "Minimum Rank-Specific Depth-of-Coverage",  value = self.configs['cutoffs']['def_val_min_rsdc'],   start = 0,   end = 1000,  step = 1),
#            "rank_level": MultiChoice(        title = "Select Taxonomic Rank to Display",         value = ["strain"],                                     options = rank_options)
#            "add_cutoff_field": Select(       title = "Select Additional Cutoff Fields:",         value = 'None',                        options = axis_options),
#            "y_axis_select": Select(          title = "Select the y-axis Value:",                 value = 'Normalized Read Count (RNR)', options = y_axis_options),
#            "add_cutoff_value": TextInput(    title = "Add an Additional Cutoff Value?"),
#            "organism_name": TextInput(       title = "Organism Name Contains ' ' :"),
#            "disease_name": TextInput(        title = "Disease Name Contains ' ' :")
            }      
            
            # Still to add:
#            "size_comparison" : Select(title="Size: Compare an Additional Measure", value='Relative Abundance', options=size_comparison_options),
#            "color_comparison" : Select(title='Color: Compare an Additional Measure', value='Score', options=color_conparison_options)
#            "patho_selector" : RadioButtonGroup(labels=["Pathogen only", "All taxonomies"], active=(0 if config['cutoffs']['def_val_patho'] else 1)),

#            rank     = Select(title='Rank', value='species', options=ranks)
#            coff_btn   = Button(label="No / Default extra cutoffs", button_type="success")
#            dl_btn  = Button(label="Export to CSV", button_type="success")

        controls_array = controls.values()

        # In documentation: the variable names within CustomJS AND all the widgets in controls{} share the names of the keys in source.data

        # IMPORTANT NOTE: ANYTHING NEW GOING INTO SOURCE MUST ALSO BE ADDED TO 'var new_data' and the '=== null' STATEMENT! DOCUMENT CAREFULLY!
        # ALSO: ENSURE THE KEY:VALUE PAIR FOR EACH WIDGET MATCHES THE NAME IN SOURCE.DATA

        # Note about MultiSelect: Just like how Read Count is on the y-axis and we see the whole axis adjust with the slider,
        # PROBABLY the MultiSelect control actually needs to be on taxa, excluding any that do not match the value(s) in multiselect.
        # I might need to, for each row, do some sort of flag and have the condition check the flag instead of the selection.

        # It's possible we'll need EVERY possible aspect of user-mediated change implemented here.
        callback = CustomJS(args = dict(source=self.source, controls=controls), code="""
            
            if (!window.full_data_save) {
                    window.full_data_save = JSON.parse(JSON.stringify(source.data));
                }

            var full_data = window.full_data_save;
            var full_data_length = full_data.taxa.length;

            var new_data = { taxa: [], taxid: [], rank_level: [], 
            raw_rd_cnt: [], tol_genome_sz: [], relative_abun: [], 
            min_score: [], score_bg: [], score_uniq: [], norm_rd_cnt: [], 
            min_norm_rd_cnt_combo: [], read_primary: [], rpkm: [], 
            linear_len: [], min_linear_cov: [], rankspec_min_depth_cov: [], min_depth_cov: [], 
            strain: [], species: [], genus: [], family: [], order: [], clade: [], phylum: [], superkingdom: [], root: [], 
            strain_rnr: [], species_rnr: [], genus_rnr: [], family_rnr: [], order_rnr: [], clade_rnr: [], phylum_rnr: [], superkingdom_rnr: [], root_rnr: [], 
            strain_rnb: [], species_rnb: [], genus_rnb: [], family_rnb: [], order_rnb: [], clade_rnb: [], phylum_rnb: [], superkingdom_rnb: [], root_rnb: [], 
            strain_ri: [], species_ri: [], genus_ri: [], family_ri: [], order_ri: [], clade_ri: [], phylum_ri: [], superkingdom_ri: [], root_ri: [], 
            patho_host: [], patho_location: [], patho_disease: [], pathogen: [], norm_value: [], NORM_REL_ABUNDANCE: []}

            console.log('')

            for (var i = 0; i < full_data_length; i++) {
                if (full_data.taxa[i] === null || full_data.taxid[i] === null || full_data.rank_level[i] === null || 
                    full_data.raw_rd_cnt === null || full_data.tol_genome_sz === null || full_data.relative_abun === null || 
                    full_data.min_score === null || full_data.score_bg === null || full_data.score_uniq === null || 
                    full_data.norm_rd_cnt === null || full_data.norm_rd_cnt_combo === null || full_data.read_primary === null || full_data.rpkm === null || 
                    full_data.linear_len === null || full_data.min_linear_cov === null || full_data.rankspec_min_depth_cov === null || full_data.min_depth_cov === null || 
                    full_data.strain === null || full_data.species === null || full_data.genus === null || full_data.family === null || full_data.order === null || full_data.clade === null || full_data.phylum === null || full_data.superkingdom === null || full_data.root === null || 
                    full_data.strain_rnr === null || full_data.species_rnr === null || full_data.genus_rnr === null || full_data.family_rnr === null || full_data.order_rnr === null || full_data.clade_rnr === null || full_data.phylum_rnr === null || full_data.superkingdom_rnr === null || full_data.root_rnr === null || 
                    full_data.strain_rnb === null || full_data.species_rnb === null || full_data.genus_rnb === null || full_data.family_rnb === null || full_data.order_rnb === null || full_data.clade_rnb === null || full_data.phylum_rnb === null || full_data.superkingdom_rnb === null || full_data.root_rnb === null || 
                    full_data.strain_ri === null || full_data.species_ri === null || full_data.genus_ri === null || full_data.family_ri === null || full_data.order_ri === null || full_data.clade_ri === null || full_data.phylum_ri === null || full_data.superkingdom_ri === null || full_data.root_ri === null ||
                    full_data.patho_host === null || full_data.patho_location === null || full_data.patho_disease === null || full_data.pathogen === null ||
                    full_data.norm_value === null || full_data.NORM_REL_ABUNDANCE === null)
                    continue;
                if (   
                    full_data.raw_rd_cnt[i] > controls.raw_rd_cnt.value &&
                    full_data.min_linear_cov[i] > controls.min_linear_cov.value &&
                    full_data.min_norm_rd_cnt_combo[i] > controls.min_norm_rd_cnt_combo.value &&
                    full_data.min_depth_cov[i] > controls.min_depth_cov.value &&
                    full_data.rankspec_min_depth_cov[i] > controls.rankspec_min_depth_cov.value &&
                    full_data.min_score[i] > controls.min_score.value
                ) { 
                    Object.keys(new_data).forEach(key => new_data[key].push(full_data[key][i]));
                }
            }

            source.data = new_data;
            source.change.emit(); 

            """)


        # Both 'controls_array' and 'callback' are referenced by routes.py, so we return both as a list to be unpacked there.
        cb_list = [controls, controls_array, callback]
        return cb_list

    ###############################################################################
    # Supporting Methods:
    ###############################################################################

    ######
    ###### INIT METHODS:
    ######

    # A method used to load a config file, which could maybe be changed in GUI settings.
    # Presently relies on hardcoding to load config_file from a JSON/Manually create config dictionary.
    def check_config(config_file = config_file):
        try:
            config_json = open(config_file, 'r').read()
            config = json.loads(config_json)
            return logconfig.bokeh_logger.debug(f"[INFO] [CONFIG] {config_file} loaded.")
        except:
            return logconfig.bokeh_logger.debug(f"[INFO] [CONFIG] {config_file} NOT loaded.")

    # A method that expects a call to Results to populate the argument with a filepath.
    def fetch_report_df(self):
        for file in os.scandir(self.filepath):
            if '.report.tsv' in Path(file).name:
                df = pd.read_csv(file, sep = '\t')
                df = df[df['LEVEL'].isin(['genus', 'species'])]
                return df

    # This method also expects a call to Results, but instead on the logfile produced in the PanGIA run.
    def fetch_log_text(self):
        for file in os.scandir(self.filepath):
            if '.pangia.log' in Path(file).name:
                with open(file, 'r') as f:
                    log = f.read()
                    return log

        ######
        ###### INIT METHODS: SIZES AND COLORS
        ######

 #   def color_fetch(self):
 #       source = self.test_source()
 #       c = [self.COLORS[8]]*len(self.report_df['NAME'].values)
 #       if source.data[colormaker] != 'None':
 #           c = [self.COLORS[int(x*110/11)] for x in self.report_df[c_col].values]
 #           return c

 #   def size_fetch(self):
 #       source = self.test_source()
 #       sz = [10]*len(self.report_df['NAME'].values)
 #       if source.data[sizemaker] != 'None':
  #          max_x = max(self.report_df[s_col].values) if max(self.report_df[s_col].values) > 0 else 0.001
  #          sz = [int(self.MIN_SIZE + x/max_x*(self.MAX_SIZE - self.MIN_SIZE)) for x in self.report_df[s_col].values]
  #          return sz

    ######
    ###### PIE CHART SUPPORTING METHODS:
    ######

    # Becaue the in-reads piechart is constructed base on information stored in the logfile 
    # (instead of the .tsv), we use this method to get counts of the reads from 'target', 'host', and 'ignored'.
    def parseLog(self):
        total_input     = 0
        total_mapped    = 0
        info            = {}
        info['Target']  = 0
        info['Host']    = 0
        info['Ignored'] = 0

        for file in os.scandir(self.filepath):
            if '.pangia.log' in Path(file).name:
                with open(file, 'r') as f:
                    for line in f:
                        if "Total number of input reads" in line:
                            reg = re.search('Total number of input reads: (\d+)', line)
                            total_input = int(reg.group(1))
                        elif "Total number of mapped reads" in line:
                            reg = re.search('Total number of mapped reads: (\d+)', line)
                            total_mapped = int(reg.group(1))
                        elif "Total number of host reads" in line:
                            reg = re.search('Total number of host reads: (\d+)', line)
                            info['Host'] = int(reg.group(1))
                        elif "Done processing SAM file" in line:
                            reg = re.search('Done processing SAM file, (\d+) alignment', line)
                            total_mapped = int(reg.group(1))
                        elif "Total number of ignored reads" in line:
                            reg = re.search('Total number of ignored reads .*: (\d+)', line)
                            info['Ignored'] = int(reg.group(1))
                            break
    
                    info['Target'] = total_mapped - info['Host'] - info['Ignored']

                    if total_input > 0:
                        info['Unmapped'] = total_input - total_mapped

                    return info

    # Called by the big pie-chart method to pull in the correct values to return to data.source;
    # This can be either from the logfile - see parseLog() - or from self.source!
    def genPieValues(self, infodict):
        tol_cnt   = 0
        ang_start = 0
        ang_size  = 6.27
        colors    = list(itemgetter(9, 3, 0, 1, 2, 5, 7, 8, 9)(Spectral11))
    
        p_name, p_str_ang, p_stp_ang, p_col, p_val, p_pct = [],[],[],[],[],[]
    
        for name in infodict:
            p_val.append(infodict[name])
            tol_cnt += infodict[name]
    
        percents = [0]
        for name in infodict:
            pct = infodict[name]/tol_cnt
            p_pct.append(pct)
            p_name.append( "%.1f%% %s"%(pct*100,name))
            percents.append( percents[-1]+(infodict[name]/tol_cnt) )
        
        p_str_ang = [p*2*pi for p in percents[:-1]]
        p_stp_ang = [p*2*pi for p in percents[1:]]
        p_col = colors[:len(p_name)]
        
        pievalue_list = [p_name, p_str_ang, p_stp_ang, p_col, p_val, p_pct]
    
        return pievalue_list


    ######
    ###### HOVERTOOL METHODS:
    ######

    # Define what the text-pop ups will say when hovering over a given value.
    # Note the one-to-one correspondance to variable names in source.data (from self.source).
    # These tools are instantiated by adding a 'tools=[]' argument to figure construction.
    def define_tooltips(self):
        hoverlist = []
        hover_graph = HoverTool(tooltips = [
                                    ("Name",                                              "@taxa"),
                                    ("TaxID",                                             "@taxid"),
                                    ("Score",                                             "@min_score{0.00}"),
                                    ("Linear Coverage",                                   "@min_linear_cov{0,0.00}"),
                                    ("Depth-of-Coverage",                                 "@min_depth_cov{0,0.00}"),
                                    ("Normalized Rank-Specific Depth-of-Coverage (RSNR)", "@rankspec_min_depth_cov{0,0.00}"),
                                    ("Relative Abundance",                                "@relative_abun{0,0.00}"),
                                    ("Raw Read Count",                                    "@raw_rd_cnt{0,0}"),
                                    ("Normalized Read Count",                             "@norm_rd_cnt{0,0.00}"),
                                    ("Normalized Rank-Specific Read Count)",              "@norm_rd_cnt_combo{0,0.00}"),
                                    ("Primary Read",                                      "@read_primary{0,0}"),
                                    ("RPKM",                                              "@rpkm{0,0.00}"),
                                    ("Score (Unique)",                                    "@score_uniq{0.00}"),
                                    ("Score (Background)",                                "@score_bg{0.00}"),
                                    ("Genome Size (bp)",                                  "@tol_genome_sz{0,0}"),
                                    ("Pathogen",                                          "@pathogen"),
                                    ])

        pie_hover = HoverTool(tooltips=[("Name",  "@taxa"),
                                        ("Reads", '@val{,} (@pct{%0.0f})' )
                                       ])


        hoverlist = [hover_graph, pie_hover]
        return hoverlist

    #######
    ####### SET UP DEFAULT CONFIGURATIONS - FIGURE SIZE & TEXT COMPONENTS
    #######

    def config_setup(self):

        # Is any of this necessary if the .json file is being loaded?!?
        config = {}
        config['cutoffs'] = {}

        config['cutoffs']['def_val_patho']         = True  #display mode: 0 -> pathogen only, 1 -> all
        config['cutoffs']['def_val_min_len']       = 50   #Minimum linear length [0-500] step=1
        config['cutoffs']['def_val_min_cov']       = 0.01  #Minimum genome coverage [0-1] step=0.01
        config['cutoffs']['def_val_max_r_raw']     = 100   #Minimum reads [0-500] step=1
        config['cutoffs']['def_val_max_r_rsnb']    = 10   #Minimum reads [0-100] step=1
        config['cutoffs']['def_val_min_score']     = 0.5   #Minimum score [0-1] step=0.1
        config['cutoffs']['def_val_min_dc']        = 10   #Minimum depth coverage MilliX (1X=1000mX) [0-10000] step=10
        config['cutoffs']['def_val_min_rsdc']      = 1     #Minimum rank specific depth coverage in MilliX[0-1000] step=0.001

        config['displays'] = {}
        config['displays']['total_width']          = 1070  #display mode: 0 -> pathogen only, 1 -> all
        config['displays']['dot_plot_width']       = 800   #Minimum linear length [0-500] step=1
        config['displays']['dot_plot_height']      = 800   #display mode: 0 -> pathogen only, 1 -> all
        config['displays']['gcov_plot_width']      = 1060  #display mode: 0 -> pathogen only, 1 -> all
        config['displays']['gcov_plot_height']     = 190   #display mode: 0 -> pathogen only, 1 -> all
        config['displays']['read_plot_width']      = 270  #display mode: 0 -> pathogen only, 1 -> all
        config['displays']['read_plot_height']     = 210   #display mode: 0 -> pathogen only, 1 -> all

        config['displays']['dashboard_pie_width']  = int(config['displays']['total_width']/2.8)
        config['displays']['dashboard_pie_height'] = int(config['displays']['dashboard_pie_width']/380*220)

        config['displays']['output_backend']       = 'canvas' #Minimum linear length [0-500] step=1

        return config
