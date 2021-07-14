import os
import ast
from datetime import datetime
import time
from flask import render_template, flash, redirect, url_for, request, g, jsonify, current_app
from flask_login import current_user, login_required
from pathlib import Path
import pandas as pd
from werkzeug.urls import url_parse
from werkzeug.utils import secure_filename
from markupsafe import Markup
from rq.job import Job

from app import db
from config import Config
from app.main.forms import AddCategoryForm, AddUserForm, AddFileTemplate, AddMetaType, StartPanGIARun, StartRealTime, ChangePanGIASettings
from app.models import User, Role, MetaTemplate, MetaType, MetaTypeToTemplates, Metadata, Category, Item, Results, ResultsToItem, Task

from app.bokeh_class import BokehObject
from bokeh.embed import components
from bokeh.resources import INLINE
from bokeh.layouts import column, row

# from app.translate import translate
# from flask_babel import _, get_locale
# from guess_language import guess_language
from app.main import bp


basedir = os.path.abspath(os.path.dirname(__file__)).replace('/app/main', '')

# Get Default PanGIA Settings
settings_file = '{}/app/static/scripts/pangia_default_params.csv'.format(basedir)
df = pd.read_csv(settings_file, sep=',', index_col='Field')
settings = df.to_dict()['Value']


@bp.before_app_request
def before_request():
    if current_user.is_authenticated:
        current_user.last_seen = datetime.utcnow()
        db.session.commit()
        # g.search_form = SearchForm()
    # g.locale = str(get_locale())


# ----------------------------------------------------------------------------------------------------------------------
# ADD EXTRA CONTEXT PROCESSORS TO TEMPLATES
# ----------------------------------------------------------------------------------------------------------------------
@bp.context_processor
def projects():
    projects = Category.query.filter_by(parent_id=0, is_archived=False, is_deleted=False)
    return {'projects': projects}


@bp.context_processor
def check_admin():
    if 'Admin' in current_user.get_roles():
        return {'is_admin': True}
    else:
        return {'is_admin': False}


# ----------------------------------------------------------------------------------------------------------------------
# DASHBOARD
# ----------------------------------------------------------------------------------------------------------------------
@bp.route('/', methods=['GET', 'POST'])
@bp.route('/index', methods=['GET', 'POST'])
@login_required
def index():
    # Beginning my attempt to replace dummy object with real one.
    # Success in retriving currently logged in user and should have access to .db objects related to current user.

    the_user = current_user
    the_user_name = (the_user.fname + ' ' + the_user.lname)
    the_user_results = Results.query.order_by(Results.results_date.desc()).limit(5)

    # Let's do the .split() on the datetime here. Keep the .html as clean as possible.

    #   pangia_runs = []
    pangia_runs_unfinished = Task.query.filter_by(complete=False).order_by(Task.task_date.desc()).all()

    #    for run in current_user.tasks:
    #        if run.get_progress == None:
    #            pangia_runs_unfinished = run.query.order_by(Task.task_date.desc()).limit(5)

    return render_template('main/dashboard.html', title='Welcome to your PanGIA Dashboard!',
                           the_user=the_user,
                           the_user_name=the_user_name,
                           the_user_results=the_user_results,
                           pangia_runs_unfinished=pangia_runs_unfinished
                           )


# ----------------------------------------------------------------------------------------------------------------------
# PROJECT PAGES
# ----------------------------------------------------------------------------------------------------------------------
@bp.route('/project/<proj_slug>', methods=['GET'])
@login_required
def project_page(proj_slug):

    # Project Not Found
    project = Category.query.filter_by(slug=proj_slug).first()
    if not project:
        return redirect(url_for('not_found_error'))

    # Get files and sort  <---- do I need to worry about pagination?
    sort_by = request.args.get('sort')
    if sort_by == 'az_up':
        files = Item.query.filter_by(category_id=project.id, is_deleted=False).order_by(Item.name.asc()).all()
    elif sort_by == 'az_down':
        files = Item.query.filter_by(category_id=project.id, is_deleted=False).order_by(Item.name.desc()).all()
    elif sort_by == 'date_up':
        files = Item.query.filter_by(category_id=project.id, is_deleted=False).order_by(Item.created_on.asc()).all()
    else:
        files = Item.query.filter_by(category_id=project.id, is_deleted=False).order_by(Item.created_on.desc()).all()

    return render_template('main/projects.html', title=project.name, project=proj_slug, files=files)


@bp.route('/project/<proj_slug>/add-file', methods=['GET', 'POST'])
@login_required
def add_file(proj_slug):

    project = Category.query.filter_by(slug=proj_slug).first()
    if not project:
        return redirect(url_for('not_found_error'))

    if request.method == 'POST':

        for f in request.files.getlist("file[]"): #request.files.items():
            # Strip out non safe characters and spaces from upload files
            filename = secure_filename(f.filename)
            filestem = Path(filename).stem
            ext = Path(filename).suffix
            # Allow uploading gzipped files
            #if ext == '.gz':
            #    ext = filename.split('.')[-2]

            # Get list of templates that include this file type
            meta_temp = MetaTemplate.query.all()

            # If extension in list it can be uploaded
            allowed_ext = []
            for x in meta_temp:
                meta_ext = str(x.file_ext).replace(' ', '').split(',')
                for y in meta_ext:
                    if y not in allowed_ext:
                        allowed_ext.append(y)

            # Upload file and put into database
            if ext in allowed_ext:
                upload_folder = '{}'.format(settings['upload_dir'])
                parent_folder = '{}{}'.format(settings['upload_dir'], proj_slug)
                upload_location = '{}{}/{}'.format(settings['upload_dir'], proj_slug, filename)
                # Make parent directory if it doesn't exist
                if not os.path.exists(upload_folder):
                    os.mkdir(upload_folder)
                if not os.path.exists(parent_folder):
                    os.mkdir(parent_folder)
                # Don't overwrite files <---------- Need to make an option so they can choose to overwrite
                if os.path.exists(upload_location):
                    flash('There is already a file with the name {} in the project {}.'.format(filename, project.name), 'warning')
                    return redirect(url_for('main.add_file', proj_slug=proj_slug))
                else:
                    f.save(upload_location)
                    add_file = Item(name=filename, item_path=upload_location, category_id=project.id, meta_template_id=0)
                    db.session.add(add_file)
                    db.session.commit()
            else:
                txt_allowed_upload = "".join(str("'" + x + "', ") for x in allowed_ext)
                flash('The only allowed files are {} files.'.format(txt_allowed_upload[:-2]), 'warning')
                return redirect(url_for('main.add_file', proj_slug=proj_slug))

            return redirect(url_for('main.add_file_meta', proj_slug=proj_slug))

    return render_template('main/add-file.html', title="Add File", project=proj_slug)


@bp.route('/project/<proj_slug>/add-file-meta', methods=['GET', 'POST'])
@login_required
def add_file_meta(proj_slug):

    project = Category.query.filter_by(slug=proj_slug).first()
    if not project:
        return redirect(url_for('not_found_error'))

    # Add metadata to file upload
    if request.method == 'POST':
        post_args = request.form
        file_id = post_args['file_id']
        template_id = post_args['template_id']
        if Item.query.filter_by(id=file_id).count() == 1:
            # Look through posted metadata and add to database <------ Need to check for required
            for x in post_args:
                if x[:4] == 'meta':
                    this_meta_id = x[5:]
                    this_meta_info = post_args[x]
                    add_meta = Metadata(item_id=file_id, metatype_id=this_meta_id, value=this_meta_info)
                    db.session.add(add_meta)
            # Update file info
            q = Item.query.filter_by(id=file_id).first()
            q.meta_template_id = template_id
            db.session.commit()

    # See if there are any files with missing meta templates if not exit
    if Item.query.filter_by(category_id=project.id, meta_template_id=0).count() < 1:
        return redirect(url_for('main.project_page', proj_slug=proj_slug))

    # Get first uploaded file that doesn't have a meta template
    file = Item.query.filter_by(category_id=project.id, meta_template_id=0).first()
    ext = Path(file.item_path).suffix
    if ext == '.gz':
        ext = str(file.item_path).split('.')[-2]

    # Get templates
    meta_templates = []
    meta_temp = MetaTemplate.query.all()
    for x in meta_temp:
        meta_ext = str(x.file_ext).replace(' ', '').split(',')
        if ext in meta_ext:
            meta_templates.append(x)

    # Get meta fields
    if len(meta_templates) == 1:
        meta_field_ids = []
        for x in meta_templates[0].metatypes:
            meta_field_ids.append(x.metatype_id)

        meta_fields = db.session.query(MetaType).filter(MetaType.id.in_(tuple(meta_field_ids))).all()

        return render_template('main/add-file_meta.html', title="Add File", project=proj_slug, file=file, meta=meta_fields,
                               meta_template=meta_templates[0].id, no_files=False)

    if len(meta_templates) > 1: # <---------------------------------------------- Need to add multiple template dropdown
        return render_template('main/add-file_meta.html', title="Add File", project=proj_slug, file=file)

    else:
        return redirect(url_for('not_found_error'))


# ----------------------------------------------------------------------------------------------------------------------
# BOKEH VISUALIZATION
# ----------------------------------------------------------------------------------------------------------------------
@bp.route('/pangia_vis/<results_id>')
def vis_button(results_id):
    # Check that there is a result with that id
    check_result = Results.query.filter_by(id=results_id).count()
    if check_result == 0:
        return redirect('main.index')
    # If not, proceed with passing plot data to components() for .html construction.
    else:
        bo = BokehObject(Results.query.get(results_id))

        # Then we pull class attributes we need to pack returned values into 'script' and 'div':
        dotplot = bo.dotplot

        pieInReads = bo.piechart_list[0]
        pieFlags = bo.piechart_list[1]
        piePatho = bo.piechart_list[2]

        datatable = bo.datatable

        # Pull the controls and CustomJS callback logic:
        cb_list = bo.control_change()
        controls_array = cb_list[1]
        callback = cb_list[2]

        # Loop through callback logic for each widget in controls array:
        for single_control in controls_array:
            single_control.js_on_change('value', callback)

        # Assemble the webpage:
        inputs_column = column(*controls_array, width=320, height=1000)
        layout_row_dotplot = row([inputs_column, dotplot])
        layout_row_piecharts = row([pieInReads, pieFlags, piePatho])
        layout_column = column([layout_row_piecharts, layout_row_dotplot, datatable])

        script, div = components(layout_column)

        # RETURN SECTION:

        # This statement can be confusing: components() RETURNS a <script> containing the data for the plots passed as arguments to the method.
        # It also returns a <div> that is the target for where the plot view is displayed.

        return render_template('items/pangia_vis.html',
                               plot_script=script, plot_div=div,
                               js_resources=INLINE.render_js(),
                               css_resources=INLINE.render_css())


# ----------------------------------------------------------------------------------------------------------------------
# ITEMS
# ----------------------------------------------------------------------------------------------------------------------
@bp.route('/item/<item_id>')
def item(item_id):
    # Check that there is an item with that id
    check_item = Item.query.filter_by(id=item_id).count()
    if check_item == 0:
        return redirect('main.index')

    # File Results
    tab = request.args.get('p')
    r_id = request.args.get('id')
    item = Item.query.get(item_id)
    if tab == 'results':

        # Go to specific results file
        if r_id and r_id != '':
            check_results = Results.query.filter_by(id=r_id).first()
            check_rtoi = ResultsToItem.query.filter_by(results_id=check_results.id, item_id=item_id).count()
            if check_results and check_rtoi > 0:
                results = Results.query.get(r_id)
                # Get File location
                for file in os.scandir(results.path):
                    if '.report.tsv' in Path(file).name:
                        tsv_file = file
                    if 'pangia.log' in Path(file).name:
                        log_file = file
                    if 'fastp.html' in Path(file).name:
                        fastp_file = Path(file).absolute()

                # Get Report TSV and format
                tsv = pd.read_csv(tsv_file, sep='\t', usecols=[
                    'LEVEL',
                    'NAME',
                    'TAXID',
                    'READ_COUNT',
                    'READ_COUNT_RNR',
                    'READ_COUNT_RSNB',
                    'LINEAR_COV',
                    'DEPTH_COV',
                    'REL_ABUNDANCE'
                ])
                tsv = tsv[tsv['LEVEL'].isin(['genus', 'species'])]
                tsv_html = Markup(tsv.to_html().replace('class="dataframe"','class="table shadow-none border-0"').replace('border="1"', '').replace('style="text-align: right;"',''))

                # Get Log file
                with open(log_file, 'r') as f:
                    log_text = f.read()

                params = ast.literal_eval(results.results_param)

                return render_template('items/item_results.html', title='Results', item=item, results=results, params=params, log_text=log_text, tsv_html=tsv_html, tsv_file=tsv_file, log_file=log_file)

        # Go to results list
        else:
            return render_template('items/item_results_list.html', title='Results List', item=item)

    # File Information
    else:
        return render_template('items/item_info.html', title='File Information', item=item)


# ----------------------------------------------------------------------------------------------------------------------
# PANGIA AND PANGIA ACTIONS
# ----------------------------------------------------------------------------------------------------------------------
@bp.route('/pangia')
def pangia():

    queue_len = Task.query.filter_by(complete=False).order_by(Task.task_date.desc()).count()
    queue = Task.query.filter_by(complete=False).order_by(Task.task_date.desc()).all()
    results = Results.query.order_by(Results.results_date.desc()).all()

    return render_template('pangia/pangia.html', title="PanGIA", queue=queue, queue_len=queue_len, results=results)


@bp.route('/pangia/processing/<task_id>')
def pangia_processing(task_id):
    t = Task.query.get(task_id)

    log_text = ''
    if t.complete == False:
        options = ast.literal_eval(t.options)
        log_file = '{}/{}.pangia.log'.format(options['tmp_dir'], options['file_stem'])
        if os.path.exists(log_file):
            with open(log_file, 'r') as f:
                log_text = f.read()

        return render_template('pangia/pangia_processing.html', task=t, options=options, log_text=log_text)
    else:
        flash('{} has finished processing.'.format(t.name), 'success')
        return redirect(url_for('main.pangia'))


@bp.route('/real_time/processing/<task_id>')
def real_time_processing(task_id):
    t = Task.query.get(task_id)

    if t.complete == False:
        options = ast.literal_eval(t.options)

        run_bokeh = True
        check_db = 0
        df = []
        try:
            # Make sure there is a result and make a Bokeh plot
            check_db = Results.query.filter_by(id=options['results_id']).count()
            result = Results.query.get(options['results_id'])
            for file in os.scandir(result.path):
                if '.report.tsv' in Path(file).name:
                    df = pd.read_csv(file, sep='\t')
        except:
            run_bokeh = False
        if check_db == 0:
            run_bokeh = False
        if len(df) == 0:
            run_bokeh = False


        if run_bokeh != False:
            bo = BokehObject(Results.query.get(options['results_id']))
            # Then we pull class attributes we need to pack returned values into 'script' and 'div':
            dotplot = bo.dotplot

            pieInReads = bo.piechart_list[0]
            pieFlags = bo.piechart_list[1]
            piePatho = bo.piechart_list[2]

            datatable = bo.datatable

            # Pull the controls and CustomJS callback logic:
            cb_list = bo.control_change()
            controls_array = cb_list[1]
            callback = cb_list[2]

            # Loop through callback logic for each widget in controls array:
            for single_control in controls_array:
                single_control.js_on_change('value', callback)

            # Assemble the webpage:
            inputs_column = column(*controls_array, width=320, height=1000)
            layout_row_dotplot = row([inputs_column, dotplot])
            layout_row_piecharts = row([pieInReads, pieFlags, piePatho])
            layout_column = column([layout_row_piecharts, layout_row_dotplot, datatable])

            script, div = components(layout_column)

            return render_template('pangia/real_time_processing.html',
                                   plot_script=script, plot_div=div, js_resources=INLINE.render_js(),
                                   css_resources=INLINE.render_css(), task=t, options=options)
        else:
            return render_template('pangia/real_time_processing.html',
                                   plot_script='', plot_div='', js_resources='',
                                   css_resources='', task=t, options=options)


    # Process Finished
    else:
        return redirect(url_for('main.pangia'))



@bp.route('/kill_job/<task_id>')
def kill_job(task_id):
    t = Task.query.get(task_id)
    t.kill_job()

    flash('{} job has been stopped.'.format(t.name), 'success')
    return redirect(url_for('main.pangia'))


@bp.route('/pangia_start', methods=['GET', 'POST'])
def pangia_start():

    global settings
    global settings_file

    for key in settings:
        if settings[key] in ['True', 'true', 'Yes', 'yes'] and settings[key] != 1:
            settings[key] = True
        if settings[key] in ['False', 'false', 'No', 'no'] and settings[key] != 0:
            settings[key] = False

    # GET ARGS to populate fields
    # Project & Item
    proj_slug = request.args.get('proj')
    item_id = request.args.get('item')
    project_list = Category.query.all()
    item = ''
    project = ''
    files = ''
    run_param = ['runname', 'description', 'seq_type', 'fastq1', 'fastq2']
    pre_param = ['preprocess', 'trim_quality', 'avg_qual_cutoff', 'min_read_len', 'n_base_cutoff', 'low_comp_filter',
                 'trim_polya', 'cut_5_prime', 'cut_3_prime']
    pangia_param = ['seed_length', 'min_align_score', 'score_method', 'min_score', 'min_read_count', 'min_read_rsnb',
                    'min_linear_len', 'min_genome_cov', 'min_depth', 'min_rs_depth', 'pathogen_discovery', 'run_annoy',
                    'run_tmark', 'run_dt']

    # Check if project is valid and get info
    if proj_slug and proj_slug != '':
        check_project = Category.query.filter_by(slug=proj_slug).count()
        if check_project > 0:
            project = Category.query.filter_by(slug=proj_slug).first()
            files = Item.query.filter_by(category_id=project.id).all()

        # Item
        if item_id and item_id != '':
            check_item = Item.query.filter_by(id=item_id, category_id=project.id).count()
            if check_item > 0:
                item = Item.query.filter_by(id=item_id, category_id=project.id).first()

        # Start form and set defaults
        form = StartPanGIARun(fastq1='fastq_{}'.format(item_id))
        form.get_fastq_choices(project.id)

        # Form submit and validation
        if form.submit.data:
            if form.validate_on_submit():
                print('valid')
                # Form is valid - redirect to the pangia action script to start the background job
                return redirect(url_for('main.action_pangia'), code=307)

        form.set_values(settings)

        all_param = run_param + pre_param + pangia_param
        errors = set()
        # Set the class for the fields based on their type and if there are errors
        for key in all_param:
            if form[key].type == 'BooleanField':
                if form[key].errors:
                    form[key].render_kw = {'class': 'form-check-input is-invalid', 'style': 'margin-left:0;'}
                    errors.add('{}<br>'.format(form[key].errors[0]))
                else:
                    form[key].render_kw = {'class': 'form-check-input', 'style': 'margin-left:0;'}
            else:
                if form[key].errors:
                    form[key].render_kw = {'class': 'form-control is-invalid'}
                    errors.add('{}<br>'.format(form[key].errors[0]))
                else:
                    form[key].render_kw = {'class': 'form-control'}

        # Flash errors to the user
        if len(errors) != 0:
            flash(Markup(''.join(errors)), 'warning')

    return render_template('pangia/run_pangia.html', title="Run PanGIA", project_list=project_list, project=project,
                           files=files, item=item, run_param=run_param, pre_param=pre_param, pangia_param=pangia_param,
                           form=form)


@bp.route('/action_pangia', methods=['POST'])
def action_pangia():

    # Format variables for pangia run
    if request.form.get('fastq1') != '' and request.form.get('fastq2') != '':
        i1 = Item.query.get(int(request.form.get('fastq1').replace('fastq_', '')))
        i2 = Item.query.get(int(request.form.get('fastq2').replace('fastq_', '')))
        fastq_loc = '{} {}'.format(i1.item_path, i2.item_path)
        items_run = [i1.id, i2.id]
    else:
        i1 = Item.query.get(int(request.form.get('fastq1').replace('fastq_', '')))
        fastq_loc = '{}'.format(i1.item_path)
        items_run = [i1.id]

    # Set up a dictionary of varibles for the pangia run
    pangia_settings = {
        'pangia_loc': settings['pangia_dir'],
        'pangia_db': settings['pangia_db'],
        'fastq_loc': fastq_loc,
        'threads': settings['threads'],
        'tmp_dir': '{}{}/results/tmp'.format(settings['upload_dir'], request.form.get('proj_slug')),
        'items_run': items_run,
        'file_stem': Path(i1.item_path).stem
    }

    run_param = ['runname', 'description', 'seq_type', 'fastq1', 'fastq2']
    pre_param = ['preprocess', 'trim_quality', 'avg_qual_cutoff', 'min_read_len', 'n_base_cutoff', 'low_comp_filter',
                 'trim_polya', 'cut_5_prime', 'cut_3_prime']
    pangia_param = ['seed_length', 'min_align_score', 'score_method', 'min_score', 'min_read_count', 'min_read_rsnb',
                    'min_linear_len', 'min_genome_cov', 'min_depth', 'min_rs_depth', 'pathogen_discovery', 'run_annoy',
                    'run_tmark', 'run_dt']

    # Add Pangia Parameters to the run settings
    for x in run_param + pre_param + pangia_param:
        pangia_settings.update({x: request.form.get(x)})

    current_user.launch_task('run_pangia', '{}'.format(request.form.get('description')), pangia_settings)
    db.session.commit()

    flash(Markup('Started a new PanGIA run for {}. See progress on <a href="{}">PanGIA results</a> page.'.format(Path(i1.item_path).name, url_for('main.pangia'))), 'success')
    sleep = time.sleep(3)
    return redirect(url_for('main.project_page', proj_slug=request.form.get('proj_slug')))


# ----------------------------------------------------------------------------------------------------------------------
# REAL TIME
# ----------------------------------------------------------------------------------------------------------------------
@bp.route('/real_time_pangia', methods=['GET', 'POST'])
def real_time():

    global settings
    global settings_file

    for key in settings:
        if settings[key] in ['True', 'true', 'Yes', 'yes'] and settings[key] != 1:
            settings[key] = True
        if settings[key] in ['False', 'false', 'No', 'no'] and settings[key] != 0:
            settings[key] = False

    # GET ARGS to populate fields
    run_param = ['runname', 'description', 'fastq_dir', 'project']
    #pre_param = ['preprocess', 'trim_quality', 'avg_qual_cutoff', 'min_read_len', 'n_base_cutoff', 'low_comp_filter',
    #             'trim_polya', 'cut_5_prime', 'cut_3_prime']
    #pangia_param = ['seed_length', 'min_align_score', 'score_method', 'min_score', 'min_read_count', 'min_read_rsnb',
    #                'min_linear_len', 'min_genome_cov', 'min_depth', 'min_rs_depth', 'pathogen_discovery', 'run_annoy',
    #                'run_tmark', 'run_dt']

    # Start form and set defaults
    form = StartRealTime()
    form.project_choices()

    # Form submit and validation
    if form.submit.data:
        if form.validate_on_submit():
            print('valid')
            # Form is valid - redirect to the pangia action script to start the background job
            return redirect(url_for('main.action_realtime'), code=307)

    #form.set_values(settings)

    all_param = run_param  # + pre_param + pangia_param
    errors = set()
    # Set the class for the fields based on their type and if there are errors
    for key in all_param:
        if form[key].type == 'BooleanField':
            if form[key].errors:
                form[key].render_kw = {'class': 'form-check-input is-invalid', 'style': 'margin-left:0;'}
                errors.add('{}<br>'.format(form[key].errors[0]))
            else:
                form[key].render_kw = {'class': 'form-check-input', 'style': 'margin-left:0;'}
        else:
            if form[key].errors:
                form[key].render_kw = {'class': 'form-control is-invalid'}
                errors.add('{}<br>'.format(form[key].errors[0]))
            else:
                form[key].render_kw = {'class': 'form-control'}

    # Flash errors to the user
    if len(errors) != 0:
        flash(Markup(''.join(errors)), 'warning')

    return render_template('pangia/real_time.html', title="Run Real Time PanGIA", run_param=run_param, form=form)


@bp.route('/action_realtime', methods=['POST'])
def action_realtime():


    # Set up a dictionary of varibles for the pangia run
    pangia_settings = {
        'pangia_loc': settings['pangia_dir'],
        'pangia_db': settings['pangia_db'],
        'threads': settings['threads'],
        'upload_dir': settings['upload_dir'],
        'fileloc': '{}{}/'.format(settings['upload_dir'], request.form.get('proj_slug'))
    }

    run_param = ['runname', 'description', 'fastq_dir', 'project']
    pre_param = ['preprocess', 'trim_quality', 'avg_qual_cutoff', 'min_read_len', 'n_base_cutoff', 'low_comp_filter',
                 'trim_polya', 'cut_5_prime', 'cut_3_prime']
    pangia_param = ['seed_length', 'min_align_score', 'score_method', 'min_score', 'min_read_count', 'min_read_rsnb',
                    'min_linear_len', 'min_genome_cov', 'min_depth', 'min_rs_depth', 'pathogen_discovery', 'run_annoy',
                    'run_tmark', 'run_dt']

    # Add Pangia Parameters to the run settings
    for x in run_param + pre_param + pangia_param:
        pangia_settings.update({x: request.form.get(x)})

    current_user.launch_task('run_real_time', '{}'.format(request.form.get('description')), pangia_settings)
    db.session.commit()

    #flash(Markup('Started a new PanGIA run for {}. See progress on <a href="{}">PanGIA results</a> page.'.format(Path(i1.item_path).name, url_for('main.pangia'))), 'success')
    sleep = time.sleep(3)
    return redirect(url_for('main.pangia'))


# ----------------------------------------------------------------------------------------------------------------------
# SETTINGS
# ----------------------------------------------------------------------------------------------------------------------
@bp.route('/settings', methods=['GET', 'POST'])
@bp.route('/settings/general', methods=['GET', 'POST'])
@login_required
def settings_page():
    if 'Admin' not in current_user.get_roles():
        return redirect(url_for('main.index'))

    global settings
    global settings_file

    for key in settings:
        if settings[key] in ['True', 'true', 'Yes', 'yes'] and settings[key] != 1:
            settings[key] = True
        if settings[key] in ['False', 'false', 'No', 'no'] and settings[key] != 0:
            settings[key] = False

    form = ChangePanGIASettings()

    # Post submit for changing pangia settings
    if form.submit.data:
        new_settings = {}
        for key in settings:
            if key in ['pangia_dir', 'pangia_db', 'upload_dir'] and form[key].data[-1] != '/':
                form[key].data = '{}/'.format(form[key].data)
            new_settings.update({key: form[key].data})
        settings = new_settings
        if form.validate_on_submit():
            new_df = pd.DataFrame({'Field': list(settings.keys()), 'Value': list(settings.values())})
            new_df.to_csv(settings_file, sep=',', index=False)
            flash('Updated settings for the PanGIA application.', 'success')

    form.set_values(settings)

    errors = set()
    # Set the class for the fields based on their type and if there are errors
    for key in settings:
        if form[key].type == 'BooleanField':
            if form[key].errors:
                form[key].render_kw = {'class': 'form-check-input is-invalid', 'style': 'margin-left:0;'}
                errors.add('{}<br>'.format(form[key].errors[0]))
            else:
                form[key].render_kw = {'class': 'form-check-input', 'style': 'margin-left:0;'}
        else:
            if form[key].errors:
                form[key].render_kw = {'class': 'form-control is-invalid'}
                errors.add('{}<br>'.format(form[key].errors[0]))
            else:
                form[key].render_kw = {'class': 'form-control'}

    # Flash errors to the user
    if len(errors) != 0:
        flash(Markup(''.join(errors)), 'warning')

    # Make a few lists to make it easier to group in html
    system_param = ['pangia_dir', 'pangia_db', 'upload_dir', 'threads']
    pre_param = ['preprocess', 'trim_quality', 'avg_qual_cutoff', 'min_read_len', 'n_base_cutoff', 'low_comp_filter',
                 'trim_polya', 'cut_5_prime', 'cut_3_prime']
    pangia_param = ['seed_length', 'min_align_score', 'score_method', 'min_score', 'min_read_count', 'min_read_rsnb',
                    'min_linear_len', 'min_genome_cov', 'min_depth', 'min_rs_depth', 'pathogen_discovery', 'run_annoy',
                    'run_tmark', 'run_dt']

    return render_template('admin/settings.html', title='General', form=form, system_param=system_param, pre_param=pre_param, pangia_param=pangia_param)


@bp.route('/settings/projects', methods=['GET', 'POST'])
@login_required
def settings_projects():
    if 'Admin' not in current_user.get_roles():
        return redirect(url_for('main.index'))

    form = AddCategoryForm(prefix='form')

    if form.submit.data:
        if form.validate_on_submit():
            # <------- Need to make these names safe
            new_cat = Category(name=form.name.data, description=form.description.data, slug=form.slug.data,
                               parent_id=form.parent_id.data, order=form.order.data)
            db.session.add(new_cat)
            db.session.commit()
            flash('Successfully added a new Project: {}'.format(form.name.data), 'success')

    categories = Category.query.all()
    return render_template('admin/settings.html', title='Projects', categories=categories, form=form)


@bp.route('/settings/templates', methods=['GET', 'POST'])
@login_required
def settings_templates():
    if 'Admin' not in current_user.get_roles():
        return redirect(url_for('main.index'))

    form = AddFileTemplate(prefix='form')
    form.metatypes.choices = form.get_metatypes()

    if form.submit.data:
        if form.validate_on_submit():
            new_template = MetaTemplate(name=form.name.data, description=form.description.data, file_ext=form.file_ext.data)
            db.session.add(new_template)
            db.session.flush()
            for x in form.metatypes.data:
                x = int(x)
                add_metatype = MetaTypeToTemplates(template_id=new_template.id, metatype_id=x, order=0, size=0)
                db.session.add(add_metatype)
            db.session.commit()

    meta_temp = MetaTemplate.query.all()
    return render_template('admin/settings.html', title='File-Templates', meta_temp=meta_temp, form=form)


@bp.route('/settings/meta-types', methods=['GET', 'POST'])
@login_required
def settings_meta_types():
    if 'Admin' not in current_user.get_roles():
        return redirect(url_for('main.index'))

    form = AddMetaType(prefix='form')

    if form.submit.data:
        if form.validate_on_submit():
            if form.value_type.data == 'choice':
                choices = form.choices.data
            else:
                choices = ''
            new_metatype = MetaType(name=form.name.data, description=form.description.data, choices=choices,
                                    value_type=form.value_type.data, is_required=form.is_required.data)
            db.session.add(new_metatype)
            db.session.commit()
            flash('Successfully added meta type - {}'.format(form.name.data), 'success')

    meta_type = MetaType.query.all()
    return render_template('admin/settings.html', title='Meta-Types', meta_type=meta_type, form=form)


# ----------------------------------------------------------------------------------------------------------------------
# USER PAGES
# ----------------------------------------------------------------------------------------------------------------------
@bp.route('/users', methods=['GET', 'POST'])
@login_required
def viewusers():
    if 'Admin' not in current_user.get_roles():
        return redirect(url_for('main.index'))

    form1 = AddUserForm(prefix='form1')
    # form2 = AddRoleForm(prefix='form2')
    # FORM 1 SUBMIT - Add new user
    if form1.submit.data:
        if form1.validate_on_submit():
            new_user = User(fname=form1.fname.data, lname=form1.lname.data, email=form1.email.data, username=form1.username.data)
            new_user.set_password(form1.password.data)
            db.session.add(new_user)
            db.session.commit()
            flash('Successfully added {} {}'.format(form1.fname.data, form1.lname.data), 'success')

    users = User.query.all()

    return render_template('admin/user.html', title='Manage Users', users=users, form1=form1)    #, form2=form2)


# Delete User ----------------------------------------------------------------------------------------------------------
@bp.route('/users/delete-user', methods=['GET', 'POST'])
@login_required
def delete_user():
    stuff = 'stuff'

