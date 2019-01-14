from flask import Blueprint, request, jsonify
# from flask_restful import Resource, Api

from app.db import db, User

bp = Blueprint('api', __name__)


# endpoint to create new user
@bp.route("/user", methods=["POST"])
def add_user():
    username = request.json['username']
    email = request.json['email']

    new_user = User(username=username, email=email)

    db.session.add(new_user)
    db.session.commit()

    return jsonify(new_user)


# endpoint to show all users
@bp.route("/user", methods=["GET"])
def get_user():
    all_users = User.query.all()

    users = {user.username: {'email': user.email} for user in all_users}
    return jsonify(users)


# endpoint to get user detail by id
@bp.route("/user/<username>", methods=["GET"])
def user_detail(username):
    user = User.query.filter_by(username=username).first_or_404()
    user = {user.username: {'email': user.email}}
    return jsonify(user)


# endpoint to update user
@bp.route("/user/<username>", methods=["PUT"])
def user_update(username):
    user = User.query.filter_by(username=username).first_or_404()
    username = request.json['username']
    email = request.json['email']

    user.email = email
    user.username = username

    db.session.commit()
    return jsonify(user)


# endpoint to delete user
@bp.route("/user/<username>", methods=["DELETE"])
def user_delete(username):
    user = User.query.filter_by(username=username).first_or_404()
    db.session.delete(user)
    db.session.commit()

    return jsonify(user)

# api = Api(bp)
#
# class TodoItem(Resource):
#     def get(self, id):
#         return {'task': 'Say "Hello, World!"'}
#
#
# api.add_resource(TodoItem, '/todos/<int:id>')


# from flask import Blueprint, send_from_directory, render_template, flash, redirect, url_for, send_file, safe_join, \
#     make_response, current_app
# from flask import jsonify
#
# import os.path as op
# from io import BytesIO
# from scipy.io import savemat, loadmat
#
# from app.models.measurement import Measurement
# from app.models.experiment import Experiment
#
# api = Blueprint('api', __name__, url_prefix='/api')
#
#
# # @api.route('/download/<string:filename>', methods=['GET', 'POST'])
# # def download(filename):
# #
# #     return send_from_directory(directory='files', filename=filename)
#
#
# @api.route('/download/matlab/<int:measurement_id>', methods=['GET', 'POST'])
# def download_matlab(measurement_id):
#     experiments = {
#         'angular': Experiment.query.get(1),
#         'spectral': Experiment.query.get(2),
#     }
#
#     measurement_data = Measurement.query.get_or_404(measurement_id)
#     file_path = op.join(current_app.root_path, 'files', measurement_data.filename)
#
#     output_data = dict()
#
#     if op.isfile(file_path):
#
#         data = loadmat(file_path, squeeze_me=True, struct_as_record=True)
#
#         if experiments['spectral'] in measurement_data.experiments \
#                 and data.keys() & {'spectral_lambda', 'spectral_I'}:
#             output_data['spectral_lambda'] = data['spectral_lambda']
#             output_data['spectral_I'] = data['spectral_I']
#
#         if experiments['angular'] in measurement_data.experiments \
#                 and data.keys() & {'angular_x', 'angular_y', 'angular_I'}:
#             output_data['angular_x'] = data['angular_x']
#             output_data['angular_y'] = data['angular_y']
#             output_data['angular_I'] = data['angular_I']
#
#     output_data['authors'] = measurement_data.authors
#
#     output_data['authors'] = measurement_data.authors
#     output_data['parameters'] = measurement_data.parameters
#     output_data['description'] = measurement_data.description
#     output_data['date_measure'] = measurement_data.date_measure.strftime('%Y-%m-%d')
#     output_data['authors'] = measurement_data.authors
#     output_data['material'] = measurement_data.material.name
#
#     buffer = BytesIO()
#     savemat(buffer, output_data)
#     buffer.seek(0)
#
#     # TODO: passing generator instead of binary data
#     return send_file(
#         buffer,
#         attachment_filename="measurement-{}-{}.mat".format(measurement_data.material.name, measurement_id),
#         as_attachment=True, conditional=True
#     )
#
#
# @api.route('/download/json/<int:measurement_id>', methods=['GET', 'POST'])
# def download_json(measurement_id):
#     experiments = {
#         'angular': Experiment.query.get(1),
#         'spectral': Experiment.query.get(2),
#     }
#
#     measurement_data = Measurement.query.get_or_404(measurement_id)
#     file_path = op.join(current_app.root_path, 'files', measurement_data.filename)
#
#     output_data = dict()
#
#     if op.isfile(file_path):
#
#         data = loadmat(file_path, squeeze_me=True, struct_as_record=True)
#
#         if experiments['spectral'] in measurement_data.experiments \
#                 and data.keys() & {'spectral_lambda', 'spectral_I'}:
#             output_data['spectral_lambda'] = data['spectral_lambda'].tolist()
#             output_data['spectral_I'] = data['spectral_I'].tolist()
#
#         if experiments['angular'] in measurement_data.experiments \
#                 and data.keys() & {'angular_x', 'angular_y', 'angular_I'}:
#             output_data['angular_x'] = data['angular_x'].tolist()
#             output_data['angular_y'] = data['angular_y'].tolist()
#             output_data['angular_I'] = data['angular_I'].tolist()
#
#     output_data['authors'] = measurement_data.authors
#
#     output_data['authors'] = measurement_data.authors
#     output_data['parameters'] = measurement_data.parameters
#     output_data['description'] = measurement_data.description
#     output_data['date_measure'] = measurement_data.date_measure.strftime('%Y-%m-%d')
#     output_data['authors'] = measurement_data.authors
#     output_data['material'] = measurement_data.material.name
#
#     output = jsonify(output_data)
#
#     output.headers['Content-Disposition'] = 'attachment;filename=measurement-{}-{}.json'.format(
#         measurement_data.material.name, measurement_id)
#
#     return output
#
#
#
# #
# # @api.route('/testt')
# # def test2():
# #     # Use BytesIO instead of StringIO here.
# #     buffer = BytesIO()
# #     buffer.write(b'jJust some letters.')
# #     # Or you can encode it to bytes.
# #     # buffer.write('Just some letters.'.encode('utf-8'))
# #     buffer.seek(0)
# #     return send_file(buffer, as_attachment=True,
# #                      attachment_filename='a_file.txt',
# #                      mimetype='text/csv')
# #
# #
# # @api.route('/test')
# # def test():
# #     csv = """"REVIEW_DATE","AUTHOR","ISBN","DISCOUNTED_PRICE"
# # "1985/01/21","Douglas Adams",0345391802,5.95
# # "1990/01/12","Douglas Hofstadter",0465026567,9.95
# # "1998/07/15","Timothy ""The Parser"" Campbell",0968411304,18.99
# # "1999/12/03","Richard Friedman",0060630353,5.95
# # "2004/10/04","Randel Helms",0879755725,4.50"""
# #     # We need to modify the response, so the first thing we
# #     # need to do is create a response out of the CSV string
# #     response = make_response(csv)
# #     # This is the key: Set the right header for the response
# #     # to be downloaded, instead of just printed on the browser
# #     response.headers["Content-Disposition"] = "attachment; filename=books.csv"
# #     return response
# #
# # from flask import stream_with_context, Response
# #
# #
# #
# # @api.route('/stream_test')
# # def stream_test():
# #     def generate():
# #         # create and return your data in small parts here
# #         for i in range(10000):
# #             yield str(i)
# #
# #     return Response(stream_with_context(generate()))
