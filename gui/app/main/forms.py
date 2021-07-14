import os
import multiprocessing
from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField, BooleanField, SubmitField, IntegerField, SelectField, SelectMultipleField, FloatField
from wtforms.validators import ValidationError, DataRequired, InputRequired, Email, EqualTo, NumberRange
from app.models import User, Role, Category, MetaType, Item
from app import db


# Custom Validators
def is_directory(form, field):
    if not os.path.exists(field.data):
        raise ValidationError('{} folder does not exist.'.format(field.label.text))


# USERS PAGE
class AddUserForm(FlaskForm):
    fname = StringField('First', validators=[DataRequired()])
    lname = StringField('Last', validators=[DataRequired()])
    username = StringField('Username', validators=[DataRequired()])
    email = StringField('Email', validators=[DataRequired(), Email()])
    password = PasswordField('Password', validators=[DataRequired()])
    submit = SubmitField('Register')

    def validate_username(self, username):
        user = User.query.filter_by(username=username.data).first()
        if user is not None:
            raise ValidationError('Please use a different username.')

    def validate_email(self, email):
        user = User.query.filter_by(email=email.data).first()
        if user is not None:
            raise ValidationError('Please use a different email address.')


# class AddRoleForm(FlaskForm):
#    rolename = StringField('Role Name', validators=[DataRequired()])
#    categories = SelectMultipleField('Select Categories', choices=get_categories())
#    can_edit = BooleanField('Can Edit')
#    can_add = BooleanField('Can Add')
#    submit = SubmitField('Add Role')
#
#    def validate_rolename(self, rolename):
#        role = Role.query.filter_by(name=rolename.data).first()
#        if role is not None:
#            raise ValidationError('There is a Role with that name.')


# SETTINGS PAGE
class ChangePanGIASettings(FlaskForm):
    pangia_dir = StringField('PanGIA Directory', validators=[DataRequired(), is_directory])
    pangia_db = StringField('PanGIA Database Location', validators=[DataRequired(), is_directory])
    upload_dir = StringField('File Upload Directory', validators=[DataRequired(), is_directory])
    threads = IntegerField('Number of Threads', validators=[NumberRange(1, multiprocessing.cpu_count(), 'Out of range')])
    preprocess = BooleanField('Run Preprocessing')
    trim_quality = IntegerField('Trim Quality Level', validators=[InputRequired()])
    avg_qual_cutoff = IntegerField('Average Quality Cutoff', validators=[InputRequired()])
    min_read_len = IntegerField('Minimum Read Length', validators=[InputRequired()])
    n_base_cutoff = IntegerField('"N" Base Cutoff', validators=[InputRequired()])
    low_comp_filter = IntegerField('Low Complexity Filter', validators=[NumberRange(0, 100, 'Out of range')])
    trim_polya = BooleanField('Trim PolyA')
    cut_5_prime = IntegerField("Cut #bp from 5'-end", validators=[InputRequired()])
    cut_3_prime = IntegerField("Cut #bp from 3'-end", validators=[InputRequired()])
    seed_length = IntegerField("Seed Length", validators=[InputRequired()])
    min_align_score = IntegerField("Minimal Aligned Score", validators=[InputRequired()])
    score_method = SelectField('Scoring Method', choices=[('standalone', 'Standalone'), ('bg', 'Background'), ('combined', 'Combined')])
    min_score = FloatField('Minimal Score', validators=[NumberRange(0, 1, 'Out of range')])
    min_read_count = IntegerField('Minimal Read Count', validators=[InputRequired()])
    min_read_rsnb = IntegerField('Minimal Read RSNB', validators=[InputRequired()])
    min_linear_len = IntegerField('Minimal Linear Length', validators=[InputRequired()])
    min_genome_cov = FloatField('Minimal Genome Coverage', validators=[InputRequired()])
    min_depth = FloatField('Minimal Depth', validators=[InputRequired()])
    min_rs_depth = FloatField('Minimal RS Depth', validators=[InputRequired()])
    pathogen_discovery = BooleanField('Pathogen Discovery')
    run_annoy = BooleanField('Run ANNOY')
    run_tmark = BooleanField('Run TMARK')
    run_dt = BooleanField('Run Decision Tree')
    # <------------------- Need decision tree options
    submit = SubmitField('Update')

    def set_values(self, settings):
        self.pangia_dir.data = settings['pangia_dir']
        self.pangia_db.data = settings['pangia_db']
        self.upload_dir.data = settings['upload_dir']
        self.threads.data = settings['threads']
        self.preprocess.data = settings['preprocess']
        self.trim_quality.data = settings['trim_quality']
        self.avg_qual_cutoff.data = settings['avg_qual_cutoff']
        self.min_read_len.data = settings['min_read_len']
        self.n_base_cutoff.data = settings['n_base_cutoff']
        self.low_comp_filter.data = settings['low_comp_filter']
        self.trim_polya.data = settings['trim_polya']
        self.cut_5_prime.data = settings['cut_5_prime']
        self.cut_3_prime.data = settings['cut_3_prime']
        self.seed_length.data = settings['seed_length']
        self.min_align_score.data = settings['min_align_score']
        self.score_method.data = settings['score_method']
        self.min_score.data = settings['min_score']
        self.min_read_count.data = settings['min_read_count']
        self.min_read_rsnb.data = settings['min_read_rsnb']
        self.min_linear_len.data = settings['min_linear_len']
        self.min_genome_cov.data = settings['min_genome_cov']
        self.min_depth.data = settings['min_depth']
        self.min_rs_depth.data = settings['min_rs_depth']
        self.pathogen_discovery.data = settings['pathogen_discovery']
        self.run_annoy.data = settings['run_annoy']
        self.run_tmark.data = settings['run_tmark']
        self.run_dt.data = settings['run_dt']


class AddCategoryForm(FlaskForm):
    name = StringField('Category Name', validators=[DataRequired('This field is required')])
    description = StringField('Description')
    slug = StringField('Slug', validators=[DataRequired('This field is required')])
    # category_type = SelectField('Category Type', choices=[('Project', 'Project')])
    parent_id = IntegerField(default=0)
    order = IntegerField(default=0)
    submit = SubmitField('Add Category')

    def validate_slug(self, slug):
        if slug.data != '':
            # ---- Make sure you strip out special characters and spaces, perferably just alphanumeric with dashes or underscores
            check = Category.query.filter_by(slug=slug.data).first()
            if check is not None:
                raise ValidationError('Categories need unique slug names.')
        else:
            raise ValidationError('Categories need unique slug names.')


class AddFileTemplate(FlaskForm):
    name = StringField('File Template Name', validators=[DataRequired('This field is required')])
    description = StringField('Description', validators=[DataRequired('This field is required')])
    file_ext = StringField('File Extensions', validators=[DataRequired('This field is required')])
    metatypes = SelectMultipleField('Select Meta Types')
    submit = SubmitField('Add Template')

    def get_metatypes(self):
        return [(str(x.id), x.name) for x in MetaType.query.all()]


class AddMetaType(FlaskForm):
    name = StringField('File Template Name', validators=[DataRequired('This field is required')])
    description = StringField('Description')
    choices = StringField('Choices')
    value_type = SelectField('Value Type', choices=[('string', 'String'), ('boolean', 'Boolean'), ('choice', 'Choice'), ('integer', 'Integer'), ('decimal', 'Decimal'), ('date', 'Date')])
    is_required = BooleanField('Is Required')
    submit = SubmitField('Add Meta Type')

    def validate_name(self, name):
        if name.data != '':
            check = MetaType.query.filter_by(name=name.data).first()
            if check is not None:
                raise ValidationError('There is already a Meta Type with that name.')
        else:
            raise ValidationError('Meta Type Name is required.')

    def validate_choices(self, choices):
        if self.value_type.data == 'choice' and str(choices.data).replace(' ', '') == '':
            raise ValidationError('You need at least one choice')


def check_float(form, field):
    try:
        field.data = float(field.data)
        if field.data > 1 or field.data < 0:
            raise ValidationError('Value must be between 0 and 1')
    except:
        raise ValidationError('Value needs to be a decimal between 0 and 1')


class StartPanGIARun(FlaskForm):
    runname = StringField('Run Name', validators=[DataRequired()])
    description = StringField('Description')
    seq_type = SelectField('Sequencer Type', choices=[('illumina', 'Illumina'), ('nanopore', 'Nanopore')], validators=[DataRequired()])
    fastq1 = SelectField('Fastq', validators=[DataRequired()])
    fastq2 = SelectField('Paired Fastq')

    preprocess = BooleanField('Run Preprocessing')
    trim_quality = IntegerField('Trim Quality Level', validators=[InputRequired()])
    avg_qual_cutoff = IntegerField('Average Quality Cutoff', validators=[InputRequired()])
    min_read_len = IntegerField('Minimum Read Length', validators=[InputRequired()])
    n_base_cutoff = IntegerField('"N" Base Cutoff', validators=[InputRequired()])
    low_comp_filter = IntegerField('Low Complexity Filter', validators=[NumberRange(0, 100, 'Out of range')])
    trim_polya = BooleanField('Trim PolyA')
    cut_5_prime = IntegerField("Cut #bp from 5'-end", validators=[InputRequired()])
    cut_3_prime = IntegerField("Cut #bp from 3'-end", validators=[InputRequired()])
    seed_length = IntegerField("Seed Length", validators=[InputRequired()])
    min_align_score = IntegerField("Minimal Aligned Score", validators=[InputRequired()])
    score_method = SelectField('Scoring Method',choices=[('standalone', 'Standalone'), ('bg', 'Background'), ('combined', 'Combined')])
    min_score = FloatField('Minimal Score', validators=[NumberRange(0, 1, 'Out of range')])
    min_read_count = IntegerField('Minimal Read Count', validators=[InputRequired()])
    min_read_rsnb = IntegerField('Minimal Read RSNB', validators=[InputRequired()])
    min_linear_len = IntegerField('Minimal Linear Length', validators=[InputRequired()])
    min_genome_cov = FloatField('Minimal Genome Coverage', validators=[InputRequired()])
    min_depth = FloatField('Minimal Depth', validators=[InputRequired()])
    min_rs_depth = FloatField('Minimal RS Depth', validators=[InputRequired()])
    pathogen_discovery = BooleanField('Pathogen Discovery')
    run_annoy = BooleanField('Run ANNOY')
    run_tmark = BooleanField('Run TMARK')
    run_dt = BooleanField('Run Decision Tree')
    submit = SubmitField('Run PanGIA')

    def validate_fastq1(self, fastq1):
        if fastq1 == '':
            raise ValidationError('You need to select a fastq to run')

    def get_fastq_choices(self, proj_id):
        q = Item.query.filter_by(category_id=proj_id).all()
        choices = [('', '-- Select Fastq --')]
        for x in q:
            choices.append(('fastq_{}'.format(x.id), x.name))
        self.fastq1.choices = choices
        self.fastq2.choices = choices

    def set_values(self, settings):
        self.preprocess.data = settings['preprocess']
        self.trim_quality.data = settings['trim_quality']
        self.avg_qual_cutoff.data = settings['avg_qual_cutoff']
        self.min_read_len.data = settings['min_read_len']
        self.n_base_cutoff.data = settings['n_base_cutoff']
        self.low_comp_filter.data = settings['low_comp_filter']
        self.trim_polya.data = settings['trim_polya']
        self.cut_5_prime.data = settings['cut_5_prime']
        self.cut_3_prime.data = settings['cut_3_prime']
        self.seed_length.data = settings['seed_length']
        self.min_align_score.data = settings['min_align_score']
        self.score_method.data = settings['score_method']
        self.min_score.data = settings['min_score']
        self.min_read_count.data = settings['min_read_count']
        self.min_read_rsnb.data = settings['min_read_rsnb']
        self.min_linear_len.data = settings['min_linear_len']
        self.min_genome_cov.data = settings['min_genome_cov']
        self.min_depth.data = settings['min_depth']
        self.min_rs_depth.data = settings['min_rs_depth']
        self.pathogen_discovery.data = settings['pathogen_discovery']
        self.run_annoy.data = settings['run_annoy']
        self.run_tmark.data = settings['run_tmark']
        self.run_dt.data = settings['run_dt']


class StartRealTime(FlaskForm):
    runname = StringField('Run Name', validators=[DataRequired()])
    description = StringField('Description')
    fastq_dir = StringField('Fastq Directory', validators=[InputRequired()])
    # output_dir = StringField('Output Directory', validators=[InputRequired()])
    project = SelectField('Project')

    #preprocess = BooleanField('Run Preprocessing')
    #trim_quality = IntegerField('Trim Quality Level', validators=[InputRequired()])
    #avg_qual_cutoff = IntegerField('Average Quality Cutoff', validators=[InputRequired()])
    #min_read_len = IntegerField('Minimum Read Length', validators=[InputRequired()])
    #n_base_cutoff = IntegerField('"N" Base Cutoff', validators=[InputRequired()])
    #low_comp_filter = IntegerField('Low Complexity Filter', validators=[NumberRange(0, 100, 'Out of range')])
    #trim_polya = BooleanField('Trim PolyA')
    #cut_5_prime = IntegerField("Cut #bp from 5'-end", validators=[InputRequired()])
    #cut_3_prime = IntegerField("Cut #bp from 3'-end", validators=[InputRequired()])
    #seed_length = IntegerField("Seed Length", validators=[InputRequired()])
    #min_align_score = IntegerField("Minimal Aligned Score", validators=[InputRequired()])
    #score_method = SelectField('Scoring Method', choices=[('standalone', 'Standalone'), ('bg', 'Background'), ('combined', 'Combined')])
    #min_score = FloatField('Minimal Score', validators=[NumberRange(0, 1, 'Out of range')])
    #min_read_count = IntegerField('Minimal Read Count', validators=[InputRequired()])
    #min_read_rsnb = IntegerField('Minimal Read RSNB', validators=[InputRequired()])
    #min_linear_len = IntegerField('Minimal Linear Length', validators=[InputRequired()])
    #min_genome_cov = FloatField('Minimal Genome Coverage', validators=[InputRequired()])
    #min_depth = FloatField('Minimal Depth', validators=[InputRequired()])
    #min_rs_depth = FloatField('Minimal RS Depth', validators=[InputRequired()])
    #pathogen_discovery = BooleanField('Pathogen Discovery')
    #run_annoy = BooleanField('Run ANNOY')
    #run_tmark = BooleanField('Run TMARK')
    #run_dt = BooleanField('Run Decision Tree')
    submit = SubmitField('Run PanGIA')

    def project_choices(self):
        q = Category.query.all()
        choices = []
        for x in q:
            choices.append((x.slug, x.name))
        self.project.choices = choices