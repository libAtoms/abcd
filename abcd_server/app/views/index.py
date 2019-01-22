from flask import Blueprint, render_template, flash, redirect, url_for, request
from flask import jsonify
import requests

from abcd_server.app.views import api

bp = Blueprint('index', __name__)

index_title = {
    'title': 'Welcome',
    'subtitle': 'ABCD database'
}

index_content = """
Lorem Ipsum is simply dummy text of the printing and typesetting industry. Lorem Ipsum has been the
standard dummy text ever since the 1500s, when an unknown printer took a galley of type and scrambled it to
make a type specimen book. It has survived not only five centuries, but also the leap into electronic
typesetting, remaining essentially unchanged. It was popularised in the 1960s with the release of Letraset
sheets containing Lorem Ipsum passages, and more recently with desktop publishing software like Aldus
PageMaker including versions of Lorem Ipsum.
"""


# Our index-page just shows a quick explanation. Check out the template
# "templates/index.html" documentation for more details.
@bp.route('/')
def index():
    return render_template("index.html",
                           title=index_title,
                           content=index_content)


@bp.route('/calculation', methods=['GET'])
def get_calculation():
    calculations = api.get_calculation().json
    return render_template("database.html",
                           title=index_title,
                           content=calculations)

