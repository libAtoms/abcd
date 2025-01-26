# Copyright (c) 2025.
# Authors: Ádám Fekete, Elliott Kasoar
# This program is distributed under the MIT License, see LICENSE.md.

from flask import Blueprint, render_template, url_for

# from flask import Blueprint, render_template, flash, redirect, url_for, request
# from flask import jsonify
# import requests

bp = Blueprint("index", __name__)


@bp.route("/")
def index():
    return render_template("index.html")


@bp.route("/login/")
def login():
    return render_template("login.html")


@bp.route("/new/")
def new():
    return render_template("new.html")


@bp.route("/graphql")
def graphql():
    return render_template("graphql.html", url=url_for("api.graphql"))
