from app import create_app, db
from app.models import User, Role,  Category, Item, MetaTemplate, MetaTypeToTemplates, MetaType, Metadata, Task, Results

app = create_app()


@app.shell_context_processor
def make_shell_context():
    return {
        'db': db,
        'User': User,
        'Role': Role,
        'Category': Category,
        'Item': Item,
        'MetaTemplate': MetaTemplate,
        'MetaTypeToTemplates': MetaTypeToTemplates,
        'MetaType': MetaType,
        'Metadata': Metadata,
        'Task': Task,
        'Results': Results
    }
