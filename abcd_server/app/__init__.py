import logging
from flask import Flask, url_for, render_template
from app.db import db

logger = logging.getLogger(__name__)



def create_app():
    app = Flask(__name__)

    app.config.from_mapping(
        SECRET_KEY='dev',
        DEBUG=True,

        # SQLALCHEMY_DATABASE_URI='sqlite:///' + os.path.join(app.root_path, 'app.db'),
        # SQLALCHEMY_DATABASE_URI='sqlite:///' + os.path.join(app.instance_path, 'app.db'),
        # SQLALCHEMY_DATABASE_URI='postgresql://docker:docker@localhost/postgres',
        # SQLALCHEMY_TRACK_MODIFICATIONS=False,

        # Configurations for PyMongo
        # MONGO_URI='mongodb://localhost:27017/abcd',
        # MONGO_URI='mongodb://fekad:qwe123@ds211613.mlab.com:11613/fekad_test'

        # Configurations for Flask_mongoengine
        MONGODB_SETTINGS={
            'db': 'atoms',
            'host': 'mongodb://localhost:27017/abcd'
        }
    )
    # app.config.from_pyfile('config.py', silent=True)
    app.config['DEBUG_TB_PANELS'] = ['flask_mongoengine.panels.MongoDebugPanel']

    # with app.app_context():
    db.init_app(app)

    # con = connect('abcd', host='mongodb://localhost/', alias='default')
    # db = con['abcd']
    # # # This will create the database file using SQLAlchemy
    # with app.app_context():
    #     db.create_all()

    # from app.views import auth
    # app.register_blueprint(auth.bp)

    from app.views import index
    app.register_blueprint(index.bp)

    from app.views import api
    app.register_blueprint(api.bp, url_prefix='/api')

    # @app.route('/user/<username>')
    # def show_user(username):
    #     user = User.query.filter_by(username=username).first_or_404()
    #     return render_template('show_user.html', user=user)

    # @app.route('/')
    # def hello_world():
    #     return 'Hello, World!'

    return app


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    from ase.io import read

    # with app.test_request_context():
    #     print(url_for('index'))
    #     print(url_for('login'))
    #     print(url_for('login', next='/'))
    #     print(url_for('profile', username='John Doe'))

    app = create_app()

    # with app.app_context():
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
# app = Flask(__name__)
# app.config.update(
#     FOO_BAR='baz',
#     FOO_SPAM='eggs',
# )
#
# foo.init_app(app)


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
# # from app.db import db
# # from app.nav import nav
# # from app.auth import basic_auth
#
# # Import database models
# # from app.models.category import Category
# # from app.models.material import Material
# # from app.models.experiment import Experiment
# # from app.models.measurement import Measurement, Experiment
# # from models.measurement import SpectralData, AngularData
#
# # Import a module / component using its blueprint handler variable (mod_auth)
# from .views.index import frontend
# # from app.views.measurements import measurements
# # from app.views.materials import materials
# # from app.views.api import api
# # from app.views.search import search
#
# # ============================================================================
# # Initialize app. Flatten config_obj to dictionary (resolve properties).
# # ============================================================================
#
# # Define the WSGI application object
# app = Flask(__name__)
#
#
# # ============================================================================
# # Initialize extensions/add-ons/plugins.
# # ============================================================================
#
# # from app.admin import admin
#
# # Configurations: The usage of Flask-Appconfig could be useful, but not a requirement.
# # app.config.from_object('config')
#
# # AppConfig(app, 'config.py')
# # Bootstrap(app)              # Install our Bootstrap extension
# # Markdown(app)
#
# # app.debug and Debug(app)    # Enable debug interface
#
#
# # db.init_app(app)            # Define the database object
# # nav.init_app(app)           # We initialize the navigation as well
# # admin.init_app(app)
# # basic_auth.init_app(app)
#
#
# # ============================================================================
# # Build the database
# # ============================================================================
#
#
# # with app.app_context():
#     # # This will create the database file using SQLAlchemy
#     # db.create_all()
#
#
# # ============================================================================
# # Setup redirects and register blueprints.
# # ============================================================================
#
#
# # app.add_url_rule('/favicon.ico', 'favicon', lambda: app.send_static_file('favicon.ico'))
#
#
# # Our application uses blueprints as well; these go well with the
# # application factory. We already imported the blueprint, now we just need
# # to register it:
# app.register_blueprint(frontend)
# # app.register_blueprint(materials)
# # app.register_blueprint(measurements)
# # app.register_blueprint(api)
# # app.register_blueprint(search)
#
# # Because we're security-conscious developers, we also hard-code disabling
# # the CDN support (this might become a default in later versions):
# # app.config['BOOTSTRAP_SERVE_LOCAL'] = True
#
#
# # Sample HTTP error handling
# @app.errorhandler(404)
# def not_found(error):
#     return render_template('404.html'), 404
#
#
# # from app.forms.search import SearchForm
# # from flask import g
#
# # @app.before_request
# # def before_request():
# #     g.search_form = SearchForm()
#
# # __all__ = [Category, Measurement, Material, db]
#
#
#
# # Migration: https://danidee10.github.io/2016/10/05/flask-by-example-5.html
