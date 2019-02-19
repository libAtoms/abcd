from flask import Blueprint, render_template, flash, redirect, url_for, request
from flask import jsonify
import requests

bp = Blueprint('database', __name__)

index_title = {
    'title': 'Welcome',
    'subtitle': 'ABCD database'
}


def get_value(row, c):
    return row[c]


@bp.route('/')
def index():
    return 200


# Our index-page just shows a quick explanation. Check out the template
# "templates/index.html" documentation for more details.
@bp.route('/<database_name>/', methods=["GET"])
def database(database_name):
    # data = Atoms.objects()
    # list(Atoms.objects.aggregate({'$unwind': '$derived.arrays_keys'}, {'$group': {'_id': '$derived.arrays_keys', 'count': {'$sum': 1}}}))

    columns = [
        {'slug': 'formula', 'name': 'Formula'},
        {'slug': 'energy', 'name': 'Energy'},
        {'slug': 'n_atoms', 'name': "# of atoms"}
    ]
    data = [
        {'id': '111', 'formula': 'Fe2O', 'energy': 1.212, 'n_atoms': 112},
        {'id': '333','formula': 'Fe2O', 'energy': 3.234, 'n_atoms': 342}
    ]

    return render_template("database/database.html",
                           database_name=database_name,
                           columns=columns,
                           data=data,
                           get_value=get_value
                           )


# Our index-page just shows a quick explanation. Check out the template
# "templates/index.html" documentation for more details.
@bp.route('/<database_name>/settings')
def settings(database_name):
    info = {
        'name': database_name,
        'description': "Lorem Ipsum is simply dummy text of the printing and typesetting industry. Lorem Ipsum has been the industry's standard dummy text ever since the 1500s, when an unknown printer took a galley of type and scrambled it to make a type specimen book. ",
        'public': True,
    }

    return render_template("database/settings.html", title=index_title, database_name=database_name, info=info)
