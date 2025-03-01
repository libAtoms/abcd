import os

from flask import Flask, render_template
from flask_nav import register_renderer

from abcd.server.app.db import db
from abcd.server.app.nav import BootstrapRenderer, DatabaseNav, nav
from abcd.server.app.views import api, database, index


def create_app(abcd_url=None):
    # Define the WSGI application object
    app = Flask(__name__)
    app.logger.info("Creating an application")

    app.config["SECRET_KEY"] = os.getenv("SECRET_KEY", os.urandom(12).hex())
    app.config["ABCD_URL"] = os.getenv("ABCD_URL", "mongodb://localhost:27017/abcd")

    # Initialize extensions/add-ons/plugins.
    nav.init_app(app)
    register_renderer(app, "BootstrapRenderer", BootstrapRenderer)
    register_renderer(app, "DatabaseNav", DatabaseNav)

    db.init_app(app)

    # Setup redirects and register blueprints.
    app.register_blueprint(index.bp)
    app.register_blueprint(database.bp, url_prefix="/db")
    app.register_blueprint(api.bp, url_prefix="/api")

    @app.errorhandler(404)
    def not_found(error):
        return render_template("404.html"), 404

    return app


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.DEBUG)

    app = create_app()
    app.run(host="0.0.0.0", debug=True)
