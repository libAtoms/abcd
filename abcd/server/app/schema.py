import graphene


class Patron(graphene.ObjectType):
    id = graphene.ID()
    name = graphene.String()
    age = graphene.Int()


class Query(graphene.ObjectType):

    patron = graphene.Field(Patron)

    def resolve_patron(self, info):
        return Patron(id=1, name="Syrus", age=27)


schema = graphene.Schema(query=Query)
query = """
    query something{
      patron {
        id
        name
        age
      }
    }
"""


def test_query():
    result = schema.execute(query)
    assert not result.errors
    assert result.data == {"patron": {"id": "1", "name": "Syrus", "age": 27}}


if __name__ == "__main__":
    result = schema.execute(query)
    print(result.data["patron"])

# import graphene
# from graphene.relay import Node
# from graphene_mongo import MongoengineConnectionField, MongoengineObjectType
#
# from app_old import models
#
#
# # #
# # class Query(graphene.ObjectType):
# #     hello = graphene.String(argument=graphene.String(default_value="stranger"))
# #
# #     def resolve_hello(self, info, argument):
# #         return 'Hello ' + argument
# #
# #
# # schema = graphene.Schema(query=Query)
#
# #
# # class Arrays(MongoengineObjectType):
# #     class Meta:
# #         model = models.Arrays
# #         interfaces = (Node,)
#
# #
# # class Info(MongoengineObjectType):
# #     class Meta:
# #         model = models.Info
# #         interfaces = (Node,)
#
# class Derived(MongoengineObjectType):
#     class Meta:
#         model = models.Derived
#
#
# class Arrays(MongoengineObjectType):
#     class Meta:
#         model = models.Arrays
#
#
# class Atoms(MongoengineObjectType):
#     class Meta:
#         model = models.Atoms
#         interfaces = (Node,)
#
#
# class Query(graphene.ObjectType):
#     node = Node.Field()
#     all_atoms = MongoengineConnectionField(Atoms)
#
#     myatoms = graphene.Field(Atoms)
#
#     def resolve_myatoms(self, info):
#         return models.Atoms.objects().first()
#
#
# #
#
# #
# # # schema = graphene.Schema(query=Query, types=[Atoms, Arrays, Info])
# schema = graphene.Schema(query=Query)
