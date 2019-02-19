import graphene
from graphene.relay import Node
from graphene_mongo import MongoengineConnectionField, MongoengineObjectType

from app import models


# #
# class Query(graphene.ObjectType):
#     hello = graphene.String(argument=graphene.String(default_value="stranger"))
#
#     def resolve_hello(self, info, argument):
#         return 'Hello ' + argument
#
#
# schema = graphene.Schema(query=Query)

#
# class Arrays(MongoengineObjectType):
#     class Meta:
#         model = models.Arrays
#         interfaces = (Node,)

#
# class Info(MongoengineObjectType):
#     class Meta:
#         model = models.Info
#         interfaces = (Node,)

class Derived(MongoengineObjectType):
    class Meta:
        model = models.Derived


class Arrays(MongoengineObjectType):
    class Meta:
        model = models.Arrays


class Atoms(MongoengineObjectType):
    class Meta:
        model = models.Atoms
        interfaces = (Node,)


class Query(graphene.ObjectType):
    node = Node.Field()
    all_atoms = MongoengineConnectionField(Atoms)

    myatoms = graphene.Field(Atoms)

    def resolve_myatoms(self, info):
        return models.Atoms.objects().first()


#

#
# # schema = graphene.Schema(query=Query, types=[Atoms, Arrays, Info])
schema = graphene.Schema(query=Query)
