import logging
from flask import Flask, render_template
from abcd_server.app.db import db

logger = logging.getLogger(__name__)


def page_not_found(e):
    return render_template('404.html'), 404


def create_app():
    app = Flask(__name__)

    app.config['DEBUG'] = True
    app.config['SECRET_KEY'] = 'sadfgsdf'

    # Configurations for Flask_mongoengine
    app.config['MONGODB_SETTINGS'] = {
        'db': 'abcd',
        'host': 'mongodb://localhost:27017/'
    }

    app.register_error_handler(404, page_not_found)

    from abcd_server.app.nav import nav, BootstrapRenderer, DatabaseNav
    from flask_nav import register_renderer
    register_renderer(app, 'BootstrapRenderer', BootstrapRenderer)
    register_renderer(app, 'DatabaseNav', DatabaseNav)
    nav.init_app(app)

    # with app_old.app_context():
    db.init_app(app)

    # con = connect('abcd', host='mongodb://localhost/', alias='default')
    # db = con['abcd']
    # # # This will create the database file using SQLAlchemy
    # with app_old.app_context():
    #     db.create_all()

    # from app_old.views import auth
    # app_old.register_blueprint(auth.bp)

    from abcd_server.app.views import index
    app.register_blueprint(index.bp)

    from abcd_server.app.views import api
    app.register_blueprint(api.bp, url_prefix='/api')

    from abcd_server.app.views import database
    app.register_blueprint(database.bp, url_prefix='/db')

    return app


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    # with app_old.test_request_context():
    #     print(url_for('index'))
    #     print(url_for('login'))
    #     print(url_for('login', next='/'))
    #     print(url_for('profile', username='John Doe'))

    app = create_app()

    # with app_old.app_context():
    #     # admin = User(username='admin', email='admin@example.com')
    #     # guest = User(username='guest', email='guest@example.com')
    #     # db.session.add(admin)
    #     # db.session.add(guest)
    #     # db.session.commit()
    #
    #     print(User.query.all())

    app.run(debug=True, host='0.0.0.0', port=5000)

    # from ase.io import write
    # from ase.io.jsonio import read_json
    # write(atoms, '-', format='json')
    #
    # traj = read('../../utils/data/bcc_bulk_54_expanded_2_high.xyz', index=slice(None))
    # for atoms in traj:
    #     # Hack to fix the representation of forces
    #     atoms.calc.results['forces'] = atoms.arrays['force']
    #
    # print(traj)
    #
    # atoms = traj[0]
    #
    # # message = json.dumps(atoms)

# # Using Extensions:

# from flask_foo import Foo
#
# foo = Foo()
#
# app_old = Flask(__name__)
# app_old.config.update(
#     FOO_BAR='baz',
#     FOO_SPAM='eggs',
# )
#
# foo.init_app(app_old)


# from flask import Flask, render_template
#
# # Import extensions
# # from flask_appconfig import AppConfig
# # from flask_bootstrap import Bootstrap
# # from flask_debug import Debug
# # from flask_markdown import Markdown
# # from flaskext.markdown import Markdown
#
#
# # from app_old.db import db
# # from app_old.nav import nav
# # from app_old.auth import basic_auth
#
# # Import database models
# # from app_old.models.category import Category
# # from app_old.models.material import Material
# # from app_old.models.experiment import Experiment
# # from app_old.models.measurement import Measurement, Experiment
# # from models.measurement import SpectralData, AngularData
#
# # Import a module / component using its blueprint handler variable (mod_auth)
# from .views.index import frontend
# # from app_old.views.measurements import measurements
# # from app_old.views.materials import materials
# # from app_old.views.api import api
# # from app_old.views.search import search
#
# # ============================================================================
# # Initialize app_old. Flatten config_obj to dictionary (resolve properties).
# # ============================================================================
#
# # Define the WSGI application object
# app_old = Flask(__name__)
#
#
# # ============================================================================
# # Initialize extensions/add-ons/plugins.
# # ============================================================================
#
# # from app_old.admin import admin
#
# # Configurations: The usage of Flask-Appconfig could be useful, but not a requirement.
# # app_old.config.from_object('config')
#
# # AppConfig(app_old, 'config.py')
# # Bootstrap(app_old)              # Install our Bootstrap extension
# # Markdown(app_old)
#
# # app_old.debug and Debug(app_old)    # Enable debug interface
#
#
# # db.init_app(app_old)            # Define the database object
# # nav.init_app(app_old)           # We initialize the navigation as well
# # admin.init_app(app_old)
# # basic_auth.init_app(app_old)
#
#
# # ============================================================================
# # Build the database
# # ============================================================================
#
#
# # with app_old.app_context():
#     # # This will create the database file using SQLAlchemy
#     # db.create_all()
#
#
# # ============================================================================
# # Setup redirects and register blueprints.
# # ============================================================================
#
#
# # app_old.add_url_rule('/favicon.ico', 'favicon', lambda: app_old.send_static_file('favicon.ico'))
#
#
# # Our application uses blueprints as well; these go well with the
# # application factory. We already imported the blueprint, now we just need
# # to register it:
# app_old.register_blueprint(frontend)
# # app_old.register_blueprint(materials)
# # app_old.register_blueprint(measurements)
# # app_old.register_blueprint(api)
# # app_old.register_blueprint(search)
#
# # Because we're security-conscious developers, we also hard-code disabling
# # the CDN support (this might become a default in later versions):
# # app_old.config['BOOTSTRAP_SERVE_LOCAL'] = True
#
#
# # Sample HTTP error handling
# @app_old.errorhandler(404)
# def not_found(error):
#     return render_template('404.html'), 404
#
#
# # from app_old.forms.search import SearchForm
# # from flask import g
#
# # @app_old.before_request
# # def before_request():
# #     g.search_form = SearchForm()
#
# # __all__ = [Category, Measurement, Material, db]
#
#
#
# # Migration: https://danidee10.github.io/2016/10/05/flask-by-example-5.html
