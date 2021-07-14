from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField, BooleanField, SubmitField
from wtforms.validators import ValidationError, DataRequired, Email, EqualTo
# from flask_babel import _, lazy_gettext as _l
from app.models import User


class LoginForm(FlaskForm):
    username = StringField('Username', validators=[DataRequired()])
    password = PasswordField('Password', validators=[DataRequired()])
    remember_me = BooleanField('Remember Me')
    submit = SubmitField('Sign In')


class FirstTimeSetup(FlaskForm):
    username = StringField('Admin Username', validators=[DataRequired()])
    password = PasswordField('Password', validators=[DataRequired()])
    fname = StringField('First Name')
    lname = StringField('Last Name')
    email = StringField('Email Address', validators=[DataRequired(), Email()])
    defaultdata = BooleanField('Install Default Setup')
    submit = SubmitField('Setup PanGIA')
