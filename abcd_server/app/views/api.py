import copy
from flask import Blueprint, request, jsonify

# from app.db import mongo



bp = Blueprint('api', __name__)


def fix_id(calc):
    calc['_id'] = str(calc['_id'])
    return calc


def make_hash(o):
    """
    Makes a hash from a dictionary, list, tuple or set to any level, that contains
    only other hashable types (including any lists, tuples, sets, and
    dictionaries).
    """

    if isinstance(o, (set, tuple, list)):

        return tuple([make_hash(e) for e in o])

    elif not isinstance(o, dict):

        return hash(o)

    new_o = copy.deepcopy(o)
    for k, v in new_o.items():
        new_o[k] = make_hash(v)

    return hash(tuple(frozenset(sorted(new_o.items()))))


class QueryString(object):
    def __init__(self, query: str):
        self.query = query

    def parse(self):
        raise NotImplementedError


# endpoint to create new user
@bp.route("/calculation", methods=["POST"])
def query_calculation():
    response = {
        'query': request.json,
        'results': []
    }
    return jsonify(response)


# endpoint to show all users
@bp.route("/calculation", methods=["GET"])
def get_calculation():

    all_calculations = [str(calculation['_id']) for calculation in mongo.db.atoms.find()]

    return jsonify(all_calculations)


# endpoint to get user detail by id
@bp.route("/calculation/<calc_id>", methods=["GET"])
def calculation_detail(calc_id):
    calculation = mongo.db.atoms.find_one_or_404(calc_id)

    return jsonify(fix_id(calculation))


# endpoint to update user
@bp.route("/calculation", methods=["PUT"])
def add_calculation():
    new_calculation = request.json
    ack = mongo.db.atoms.insert_one(new_calculation)

    # return make_response('', 403)
    return jsonify(str(ack.inserted_id))


# endpoint to delete user
@bp.route("/calculation/<calc_id>", methods=["DELETE"])
def delete_calculation(calc_id):
    calculation = mongo.db.atoms.find_one_or_404(calc_id)
    _ = mongo.db.atoms.delete_one(calculation)

    return jsonify(fix_id(calculation))
    # return make_response('The item has been deleted!', 200)


@bp.route("/search/", methods=['GET'])
def search_calculation():
    data = [str(calc['_id']) for calc in mongo.db.atoms.find()]
    return jsonify(data)

# # @api.route('/stream_test')
# # def stream_test():
# #     def generate():
# #         # create and return your data in small parts here
# #         for i in range(10000):
# #             yield str(i)
# #
# #     return Response(stream_with_context(generate()))
