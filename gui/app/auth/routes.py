from flask import render_template, redirect, url_for, flash, request, current_app
from werkzeug.urls import url_parse
from flask_login import login_user, logout_user, current_user
# from flask_babel import _
import pandas as pd
from sqlalchemy import DDL
from app import db
from app.auth import bp
from app.auth.forms import LoginForm, FirstTimeSetup
from app.models import User, Notification, Role, MetaTemplate, MetaType, MetaTypeToTemplates, Metadata, Category, Item, Results, ResultsToItem, Task
from app.auth.email import send_password_reset_email


# ----------------------------------------------------------------------------------------------------------------------
# SITE LOGIN
# ----------------------------------------------------------------------------------------------------------------------
@bp.route('/login', methods=['GET', 'POST'])
def login():
    # Check to see if there are any users
    # If no users than it could be a first time setup
    num_users = User.query.count()
    if num_users == 0:
        return redirect(url_for('auth.first_time_setup'))

    # Login Form
    else:

        if current_user.is_authenticated:
            return redirect(url_for('main.index'))
        form = LoginForm()
        if form.validate_on_submit():
            user = User.query.filter_by(username=form.username.data).first()
            if user is None or not user.check_password(form.password.data):
                flash('Invalid username or password')
                return redirect(url_for('auth.login'))
            login_user(user, remember=form.remember_me.data)
            next_page = request.args.get('next')
            if not next_page or url_parse(next_page).netloc != '':
                next_page = url_for('main.index')
            return redirect(next_page)
        return render_template('auth/login.html', title='Sign In', form=form)


@bp.route('/logout')
def logout():
    logout_user()
    return redirect(url_for('auth.login'))


@bp.route('/register', methods=['GET', 'POST'])
def register():
    return True


@bp.route('/reset_password_request', methods=['GET', 'POST'])
def reset_password_request():
    return True


@bp.route('/reset_password/<token>', methods=['GET', 'POST'])
def reset_password(token):
    return True

# ----------------------------------------------------------------------------------------------------------------------
# DEBUG SITE RESET - Comment out after development
# ----------------------------------------------------------------------------------------------------------------------
@bp.route('/reset_site')
def reset_site():

    # Delete all data
    User.query.delete()
    Notification.query.delete()
    Role.query.delete()
    Category.query.delete()
    Item.query.delete()
    Results.query.delete()
    ResultsToItem.query.delete()
    Metadata.query.delete()
    MetaTemplate.query.delete()
    MetaType.query.delete()
    MetaTypeToTemplates.query.delete()
    Task.query.delete()
    db.session.commit()

    # Reset auto increment
    DDL("ALTER TABLE {} AUTO_INCREMENT = 1".format(User.__table__))
    DDL("ALTER TABLE {} AUTO_INCREMENT = 1".format(Notification.__table__))
    DDL("ALTER TABLE {} AUTO_INCREMENT = 1".format(Role.__table__))
    DDL("ALTER TABLE {} AUTO_INCREMENT = 1".format(Category.__table__))
    DDL("ALTER TABLE {} AUTO_INCREMENT = 1".format(Item.__table__))
    DDL("ALTER TABLE {} AUTO_INCREMENT = 1".format(Results.__table__))
    DDL("ALTER TABLE {} AUTO_INCREMENT = 1".format(ResultsToItem.__table__))
    DDL("ALTER TABLE {} AUTO_INCREMENT = 1".format(Metadata.__table__))
    DDL("ALTER TABLE {} AUTO_INCREMENT = 1".format(MetaTemplate.__table__))
    DDL("ALTER TABLE {} AUTO_INCREMENT = 1".format(MetaType.__table__))
    DDL("ALTER TABLE {} AUTO_INCREMENT = 1".format(MetaTypeToTemplates.__table__))
    DDL("ALTER TABLE {} AUTO_INCREMENT = 1".format(Task.__table__))

    return redirect(url_for('main.index'))


# ----------------------------------------------------------------------------------------------------------------------
# FIRST TIME GUI SETUP
# ----------------------------------------------------------------------------------------------------------------------
@bp.route('/setup', methods=['GET', 'POST'])
def first_time_setup():
    # Check that there hasn't already been an admin made
    # If no admin then create one
    if (Role.query.filter_by(role_type='Admin').count() == 0) and (User.query.count() == 0):
        form = FirstTimeSetup()

        if form.validate_on_submit():

            # Add the administrative user
            new_user = User(fname=form.fname.data, lname=form.lname.data, email=form.email.data, username=form.username.data)
            new_user.set_password(form.password.data)
            db.session.add(new_user)
            db.session.flush()

            # Add admin role
            admin_role = Role(role_type='Admin', user_id=new_user.id)
            db.session.add(admin_role)

            # Read from the default setup file and add to database
            df = pd.read_csv('{}/app/static/scripts/pangia_default_db.csv'.format(current_app.config['BASEDIR']), sep=',')

            # Iterate through Templates
            meta_templates = df[df['Type'] == 'MetaTemplate']
            for key, row in meta_templates.iterrows():

                # Add meta templates
                new_metatemplate = MetaTemplate(name=row['Name'], description=row['Description'], file_ext=row['Value'])
                db.session.add(new_metatemplate)
                db.session.flush()

                # Iterate through Meta and add meta types
                this_metatypes = df[df['Associations'] == row['Name']]
                for mkey, mrow in this_metatypes.iterrows():

                    # Get variables
                    is_bool = False
                    is_choice = False
                    is_required = False
                    choices = ''
                    if mrow['Value'] == 'choice':
                        choices = mrow['Options']
                        is_choice = 1
                    elif mrow['Value'] == 'boolean':
                        is_bool = 1
                    if row['Is_Required'] == 'True':
                        is_required = True

                    # Make new metatype
                    new_metatype = MetaType(name=mrow['Name'], description=mrow['Description'], value_type=mrow['Value'],
                                            choices=choices, is_bool=is_bool, is_choice=is_choice,
                                            is_required=is_required)
                    db.session.add(new_metatype)
                    db.session.flush()

                    # Make metatype to template association
                    new_metatotemplate = MetaTypeToTemplates(template_id=new_metatemplate.id, metatype_id=new_metatype.id, order=0)
                    db.session.add(new_metatotemplate)

            db.session.commit()
            return render_template('auth/setup.html', title='Setup Page', complete=True)

        return render_template('auth/setup.html', tilte='Setup Page', form=form, complete=False)

    # Already an admin, redirect to home page
    else:
        return redirect(url_for('main.index'))
