from flask import jsonify
from flask import Blueprint, render_template, flash, redirect, url_for, request
import requests

from mongoengine.context_managers import switch_collection

# from abcd_server.db.backends.mongoengine import Atoms


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

    db_info = {
        'name': database_name,
        'description': 'Vivamus sagittis lacus vel augue laoreet rutrum faucibus dolor auctor. Duis mollis, est non commodo luctus.',
        'columns': [
            {'slug': 'formula', 'name': 'Formula'},
            {'slug': 'energy', 'name': 'Energy'},
            {'slug': 'n_atoms', 'name': "# of atoms"}
        ],
    }

    columns = [
        {'slug': 'formula', 'name': 'Formula'},
        {'slug': 'energy', 'name': 'Energy'},
        {'slug': 'n_atoms', 'name': "# of atoms"}
    ]

    data = [
        {'id': '111', 'formula': 'Fe2O', 'energy': 1.212, 'n_atoms': 112},
        {'id': '333', 'formula': 'Fe2O', 'energy': 3.234, 'n_atoms': 342}
    ]

    from abcd_server.app.models import Atoms
    page = request.args.get('page', 1, type=int)

    with switch_collection(Atoms, database_name) as Atoms:
        # atoms = Atoms.objects.get_or_404()
        atoms = Atoms.objects
        paginated_atoms = atoms.paginate(page, per_page=10)
        # paginated_atoms = atoms.paginate_field('tags', page, per_page=10)

        # print(Atoms.objects.paginate(page=page, per_page=10))

    # list(Atoms.objects.aggregate({'$unwind': '$derived.arrays_keys'}, {'$group': {'_id': '$derived.arrays_keys', 'count': {'$sum': 1}}}))

    return render_template("database/database.html",
                           db_info=db_info,
                           columns=columns,
                           data=data,
                           atoms=paginated_atoms)


# Our index-page just shows a quick explanation. Check out the template
# "templates/index.html" documentation for more details.
@bp.route('/<database_name>/settings')
def settings(database_name):
    info = {
        'name': database_name,
        'description': "Lorem Ipsum is simply dummy text of the printing and typesetting industry. Lorem Ipsum has been the industry's standard dummy text ever since the 1500s, when an unknown printer took a galley of type and scrambled it to make a type specimen book. ",
        'public': True,
    }

    return render_template("database/settings.html", database_name=database_name, info=info)
