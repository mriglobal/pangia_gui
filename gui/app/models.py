from datetime import datetime
from hashlib import md5
import json
from time import time
from flask import current_app
from flask_login import UserMixin
from werkzeug.security import generate_password_hash, check_password_hash
import jwt
import redis
import rq
from rq.command import send_stop_job_command
from app import db, login
#from app.search import add_to_index, remove_from_index, query_index


'''
class SearchableMixin(object):
    @classmethod
    def search(cls, expression, page, per_page):
        ids, total = query_index(cls.__tablename__, expression, page, per_page)
        if total == 0:
            return cls.query.filter_by(id=0), 0
        when = []
        for i in range(len(ids)):
            when.append((ids[i], i))
        return cls.query.filter(cls.id.in_(ids)).order_by(
            db.case(when, value=cls.id)), total

    @classmethod
    def before_commit(cls, session):
        session._changes = {
            'add': list(session.new),
            'update': list(session.dirty),
            'delete': list(session.deleted)
        }

    @classmethod
    def after_commit(cls, session):
        for obj in session._changes['add']:
            if isinstance(obj, SearchableMixin):
                add_to_index(obj.__tablename__, obj)
        for obj in session._changes['update']:
            if isinstance(obj, SearchableMixin):
                add_to_index(obj.__tablename__, obj)
        for obj in session._changes['delete']:
            if isinstance(obj, SearchableMixin):
                remove_from_index(obj.__tablename__, obj)
        session._changes = None

    @classmethod
    def reindex(cls):
        for obj in cls.query:
            add_to_index(cls.__tablename__, obj)


db.event.listen(db.session, 'before_commit', SearchableMixin.before_commit)
db.event.listen(db.session, 'after_commit', SearchableMixin.after_commit)
'''


class User(db.Model, UserMixin):
    __tablename__ = 'user'
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(64), index=True, unique=True)
    fname = db.Column(db.String(20))
    lname = db.Column(db.String(20))
    email = db.Column(db.String(120), index=True, unique=True)
    password_hash = db.Column(db.String(128))
    last_seen = db.Column(db.DateTime, default=datetime.utcnow)
    roles = db.relationship('Role', backref='user', lazy='dynamic')
    notifications = db.relationship('Notification', backref='user', lazy='dynamic')
    tasks = db.relationship('Task', backref='user', lazy='dynamic')
    results = db.relationship('Results', backref='user', lazy='dynamic')

    def __repr__(self):
        return '<User {}>'.format(self.username)

    def set_password(self, password):
        self.password_hash = generate_password_hash(password)

    def check_password(self, password):
        return check_password_hash(self.password_hash, password)

    def get_reset_password_token(self, expires_in=600):
        return jwt.encode(
            {'reset_password': self.id, 'exp': time() + expires_in},
            current_app.config['SECRET_KEY'],
            algorithm='HS256').decode('utf-8')

    def get_roles(self):
        roles = []
        for x in Role.query.filter_by(user_id=self.id).all():
            roles.append(x.role_type)
        return roles

    @staticmethod
    def verify_reset_password_token(token):
        try:
            id = jwt.decode(token, current_app.config['SECRET_KEY'],
                            algorithms=['HS256'])['reset_password']
        except:
            return
        return User.query.get(id)

    def add_notification(self, name, data):
        self.notifications.filter_by(name=name).delete()
        n = Notification(name=name, payload_json=json.dumps(data), user=self)
        db.session.add(n)
        return n

    def launch_task(self, name, description, args):
        rq_job = current_app.task_queue.enqueue('app.tasks.' + name, self.id, args, job_timeout=3600)
        #try:
        task = Task(id=rq_job.get_id(), task_type=name, name=args['runname'],
                    description=description, options=str(args), user=self)
        #except KeyError:
        #    task = Task(id=rq_job.get_id(), task_type=name, description=description, user=self)
        db.session.add(task)
        return task

    def get_tasks_in_progress(self):
        return Task.query.filter_by(user=self, complete=False).all()

    def get_task_in_progress(self, name):
        return Task.query.filter_by(name=name, user=self, complete=False).first()


@login.user_loader
def load_user(id):
    return User.query.get(int(id))


class Notification(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(128), index=True)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'))
    timestamp = db.Column(db.Float, index=True, default=time)
    payload_json = db.Column(db.Text)

    def get_data(self):
        return json.loads(str(self.payload_json))


# ----------------------------------------------------------------------------------------------------------------------
# USER PERMISIONS
# ----------------------------------------------------------------------------------------------------------------------
# I ended up nerfing Roles a lot and removing most of its functionality
class Role(db.Model):
    __tablename__ = 'roles'
    id = db.Column(db.Integer(), primary_key=True)
    role_type = db.Column(db.String(128), unique=True)
    user_id = db.Column(db.Integer(), db.ForeignKey('user.id', ondelete='CASCADE'))
    # special = db.Column(db.String(255))
    # can_edit = db.Column(db.Boolean, default=False)
    # can_add = db.Column(db.Boolean, default=False)
    # Relationships
    # category = db.relationship('CategoryRoles', foreign_keys='CategoryRoles.role_id')
    # users = db.relationship('UserRoles', foreign_keys='user_roles.role_id')


# class UserRoles(db):
#    __tablename__ = 'user_roles'
#    id = db.Column(db.Integer(), primary_key=True)
#    user_id = db.Column(db.Integer(), db.ForeignKey('user.id', ondelete='CASCADE'))
#    role_id = db.Column(db.Integer(), db.ForeignKey('roles.id', ondelete='CASCADE'))
#    username = db.relationship('User', foreign_keys="UserRoles.user_id")


# class CategoryRoles(db):
#    __tablename__ = 'category_roles'
#    id = db.Column(db.Integer(), primary_key=True)
#    role_id = db.Column(db.Integer(), db.ForeignKey('roles.id', ondelete='CASCADE'))
#    category_id = db.Column(db.Integer(), db.ForeignKey('category.id', ondelete='CASCADE'))
#    name = db.relationship('Category', foreign_keys='CategoryRoles.category_id')


# ----------------------------------------------------------------------------------------------------------------------
# PROJECTS AND PAGES
# ----------------------------------------------------------------------------------------------------------------------
class Category(db.Model):
    __tablename__ = 'category'
    id = db.Column(db.Integer(), primary_key=True)
    name = db.Column(db.String(128), index=True)
    description = db.Column(db.Text)
    slug = db.Column(db.String(20), index=True)
    # category_type = db.Column(db.String(128))
    parent_id = db.Column(db.Integer(), db.ForeignKey('category.id', name='parent_id'))
    root_id = db.Column(db.Integer(), index=True, default=0)
    order = db.Column(db.Integer(), default=0)
    created_on = db.Column(db.DateTime, default=datetime.utcnow)
    is_archived = db.Column(db.Boolean, default=False)
    is_deleted = db.Column(db.Boolean, default=False)
    # items = db.relationship('Item', foreign_keys="[Item.category_id]")
    # meta = db.relationship('CategoryMeta')
    items = db.relationship('Item', backref='category', lazy='dynamic')
    children = db.relationship('Category', foreign_keys="Category.parent_id")

    def __repr__(self):
        return '<Category {}>'.format(self.name)


# class CategoryMeta(db):
#    __tablename__ = 'category_meta'
#    id = db.Column(db.Integer(), primary_key=True)
#    name = db.Column(db.String(255), unique=True, index=True)
#    value = db.Column(db.Text)
#    category_id = db.Column(db.Integer, db.ForeignKey('category.id'))
#
#    def __repr__(self):
#        return '<CategoryMeta {}>'.format(self.name)


# ----------------------------------------------------------------------------------------------------------------------
# FILES
# ----------------------------------------------------------------------------------------------------------------------
class Item(db.Model):
    __tablename__ = 'item'
    id = db.Column(db.Integer(), primary_key=True)
    name = db.Column(db.String(128), index=True)
    item_path = db.Column(db.String(128))
    description = db.Column(db.Text)
    category_id = db.Column(db.Integer(), db.ForeignKey('category.id'))
    meta_template_id = db.Column(db.Integer, db.ForeignKey('meta_template.id'))
    created_on = db.Column(db.DateTime, index=True, default=datetime.utcnow)
    is_archived = db.Column(db.Boolean, default=False)
    is_deleted = db.Column(db.Boolean, default=False)
    # Relationships
    meta_template = db.relationship('MetaTemplate', foreign_keys=[meta_template_id])
    metas = db.relationship('Metadata', foreign_keys="[Metadata.item_id]")
    results = db.relationship('Results', secondary='results_to_item')

    def __repr__(self):
        return '<Item {}>'.format(self.name)


# ----------------------------------------------------------------------------------------------------------------------
# RESULTS
# ----------------------------------------------------------------------------------------------------------------------
class Results(db.Model):
    __tablename__ = 'results'
    id = db.Column(db.Integer(), primary_key=True)
    name = db.Column(db.String(128), index=True)
    description = db.Column(db.String(128))
    path = db.Column(db.String(128))
    results_type = db.Column(db.String(128))
    results_date = db.Column(db.DateTime, index=True, default=datetime.utcnow)
    results_param = db.Column(db.Text)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'))
    items = db.relationship('Item', secondary='results_to_item')


class ResultsToItem(db.Model):
    __tablename__ = 'results_to_item'
    results_id = db.Column(db.Integer, db.ForeignKey('results.id'), primary_key=True)
    item_id = db.Column(db.Integer, db.ForeignKey('item.id'), primary_key=True)


# ----------------------------------------------------------------------------------------------------------------------
# FILE METADATA
# ----------------------------------------------------------------------------------------------------------------------
class MetaTemplate(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(128))
    description = db.Column(db.Text)
    file_ext = db.Column(db.String(128))
    # Relationships
    metatypes = db.relationship('MetaTypeToTemplates', foreign_keys="[MetaTypeToTemplates.template_id]")

    def __repr__(self):
        return '<MetaTemplate {}>'.format(self.name)


class MetaTypeToTemplates(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    template_id = db.Column(db.Integer, db.ForeignKey('meta_template.id', ondelete='CASCADE'))
    metatype_id = db.Column(db.Integer, db.ForeignKey('meta_type.id', ondelete='CASCADE'))
    order = db.Column(db.Integer)
    size = db.Column(db.Integer)  # I don't use this. Was going to be the size of the input field

    def get_name(self):
        query = db.session.query(MetaType.name).filter_by(id=self.metatype_id).first()
        return query.name
    def get_id(self):
        query = db.session.query(MetaType.id).filter_by(id=self.metatype_id).first()
        return query.id


class MetaType(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(128))
    description = db.Column(db.Text)
    value_type = db.Column(db.String(128), default="string")
    choices = db.Column(db.Text)
    is_bool = db.Column(db.Boolean, default=0)
    is_choice = db.Column(db.Boolean, default=0)
    is_required = db.Column(db.Boolean, default=0)

    def __repr__(self):
        return '<MetaType {}>'.format(self.name)

    # def get_values(self):
    #    values = []
    #    query = db.session.query(Metadata.value).filter(Metadata.metatype_id==self.id).all()
    #    return Counter(query)


class Metadata(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    item_id = db.Column(db.Integer, db.ForeignKey('item.id'), index=True)
    metatype_id = db.Column(db.Integer, db.ForeignKey('meta_type.id'))
    value = db.Column(db.Text)
    # Relationships
    metatype = db.relationship('MetaType', foreign_keys=[metatype_id])

    def __repr__(self):
        return '<Metadata {}>'.format(self.value)


# ----------------------------------------------------------------------------------------------------------------------
# BACKGROUND TASKS
# ----------------------------------------------------------------------------------------------------------------------
class Task(db.Model):
    id = db.Column(db.String(36), primary_key=True)
    name = db.Column(db.String(128), index=True)
    task_type = db.Column(db.String(128), index=True)
    description = db.Column(db.String(128))
    options = db.Column(db.String(1024))
    task_date = db.Column(db.DateTime, index=True, default=datetime.utcnow)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'))
    complete = db.Column(db.Boolean, default=False)

    def get_rq_job(self):
        try:
            rq_job = rq.job.Job.fetch(self.id, connection=current_app.redis)
        except (redis.exceptions.RedisError, rq.exceptions.NoSuchJobError):
            return None
        return rq_job

    def get_progress(self):
        job = self.get_rq_job()
        return job.meta.get('progress', 0) if job is not None else 100

    def kill_job(self):
        try:
            send_stop_job_command(connection=current_app.redis, job_id=self.id)
        except:
            not_running = True
        self.complete = True
        db.session.commit()