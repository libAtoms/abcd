from flask import Blueprint, render_template, request
from flask import Response

bp = Blueprint("database", __name__)


@bp.route("/")
def index():
    return Response(status=200)


# Our index-page just shows a quick explanation. Check out the template
# "templates/index.html" documentation for more details.
@bp.route("/<database_name>/", methods=["GET"])
def database(database_name):
    # data = Atoms.objects()
    # list(Atoms.objects.aggregate({'$unwind': '$derived.arrays_keys'}, {'$group': {'_id': '$derived.arrays_keys', 'count': {'$sum': 1}}}))

    if request.method == "POST":
        print("POST")

    info = {
        "name": database_name,
        "description": "Vivamus sagittis lacus vel augue laoreet rutrum faucibus dolor auctor. Duis mollis, est non commodo luctus.",
        "columns": [
            {"slug": "formula", "name": "Formula"},
            {"slug": "energy", "name": "Energy"},
            {"slug": "derived.n_atoms", "name": "# of atoms"},
        ],
    }

    paginated_atoms = []

    # page = request.args.get('page', 1, type=int)
    # paginated_atoms = atoms.paginate(page, per_page=10)

    return render_template("database/database.html", info=info, atoms=paginated_atoms)


# Our index-page just shows a quick explanation. Check out the template
# "templates/index.html" documentation for more details.
@bp.route("/<database_name>/settings")
def settings(database_name):
    info = {
        "name": database_name,
        "description": "Vivamus sagittis lacus vel augue laoreet rutrum faucibus dolor auctor. Duis mollis, est non commodo luctus.",
        "columns": [
            {"slug": "formula", "name": "Formula"},
            {"slug": "energy", "name": "Energy"},
            {"slug": "derived.n_atoms", "name": "# of atoms"},
        ],
        "public": True,
    }
    return render_template(
        "database/settings.html", database_name=database_name, info=info
    )
