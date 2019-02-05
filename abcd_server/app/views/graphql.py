from flask import Blueprint, request, jsonify
from flask_graphql import GraphQLView

from app.schema import schema

bp = Blueprint('graphql', __name__)

bp.add_url_rule('/', view_func=GraphQLView.as_view('graphql', schema=schema, graphiql=True))

#
# @bp.route('/')
# def index():
#     return jsonify('hello')

# from database import init_db
# from flask import Flask
# from flask_graphql import GraphQLView
# from app.schema import schema
#
# default_query = '''
# {
#   allEmployees {
#     edges {
#       node {
#         id,
#         name,
#         department {
#           id,
#           name
#         },
#         roles {
#           edges {
#             node {
#               id,
#               name
#             }
#           }
#         },
#         leader {
#           id,
#           name
#         }
#         tasks {
#           edges {
#             node {
#               name,
#               deadline
#             }
#           }
#         }
#       }
#     }
#   }
# }'''.strip()
#
