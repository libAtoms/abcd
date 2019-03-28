from flask import Blueprint, render_template, request

from mongoengine.context_managers import switch_collection

# from server.db.backends.mongoengine import Atoms


bp = Blueprint('database', __name__)


@bp.route('/')
def index():
    return 200


# Our index-page just shows a quick explanation. Check out the template
# "templates/index.html" documentation for more details.
@bp.route('/<database_name>/', methods=['GET'])
def database(database_name):
    # data = Atoms.objects()
    # list(Atoms.objects.aggregate({'$unwind': '$derived.arrays_keys'}, {'$group': {'_id': '$derived.arrays_keys', 'count': {'$sum': 1}}}))
    if request.method == 'POST':
        print('POST')

    print(request.method)

    info = {
        'name': database_name,
        'description': 'Vivamus sagittis lacus vel augue laoreet rutrum faucibus dolor auctor. Duis mollis, est non commodo luctus.',
        'columns': [
            {'slug': 'formula', 'name': 'Formula'},
            {'slug': 'energy', 'name': 'Energy'},
            {'slug': 'derived.n_atoms', 'name': "# of atoms"}
        ],
    }

    from abcd.server.app.models import Atoms
    page = request.args.get('page', 1, type=int)

    with switch_collection(Atoms, database_name) as Atoms:
        # atoms = Atoms.objects.get_or_404()
        atoms = Atoms.objects
        paginated_atoms = atoms.paginate(page, per_page=10)
        # paginated_atoms = atoms.paginate_field('tags', page, per_page=10)

    # list(Atoms.objects.aggregate({'$unwind': '$derived.arrays_keys'}, {'$group': {'_id': '$derived.arrays_keys', 'count': {'$sum': 1}}}))

    return render_template("database/database.html",
                           info=info,
                           atoms=paginated_atoms)


# Our index-page just shows a quick explanation. Check out the template
# "templates/index.html" documentation for more details.
@bp.route('/<database_name>/settings')
def settings(database_name):

    info = {
        'name': database_name,
        'description': 'Vivamus sagittis lacus vel augue laoreet rutrum faucibus dolor auctor. Duis mollis, est non commodo luctus.',
        'columns': [
            {'slug': 'formula', 'name': 'Formula'},
            {'slug': 'energy', 'name': 'Energy'},
            {'slug': 'derived.n_atoms', 'name': "# of atoms"}
        ],
        'public': True
    }
    return render_template("database/settings.html", database_name=database_name, info=info)
