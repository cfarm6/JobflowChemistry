from typing import Dict, Iterable, List
from icecream import ic
from maggma.stores import MongoStore
from maggma.core import Builder, Store
import logging
#
workflow_store = MongoStore(
    database="fireworks",
    collection_name="workflows"
)
fireworks_store = MongoStore(
        database="fireworks",
        collection_name="fireworks"
        )
docs_store = MongoStore(database="data", collection_name="docs_store")
property_store = MongoStore(
    database="data",
    collection_name="properties"
)
molecule_property_store = MongoStore(
    database="molecule",
    collection_name="property"
)
# Builder
class MolecularProperties(Builder):
    fields: List = ["global", "atomic", "bonds"]
    def __init__(
        self,
        workflow_store: Store,
        fireworks_store: Store,
        docs_store: Store,
        property_store: Store,
        target_store: Store,
    ):
        self.workflow_store = workflow_store
        self.fireworks_store = fireworks_store
        self.docs_store = docs_store
        self.property_store = property_store
        self.target_store = target_store
        super().__init__(sources=[workflow_store, fireworks_store, docs_store, property_store], targets = target_store)

    def get_items(self) -> Iterable:
        workflows = self.workflow_store.query()
        for workflow in workflows:
            cid = int(workflow['name'])
            fireworks = fireworks_store.query({"fw_id": {"$in": workflow["nodes"]}})
            for firework in fireworks:
                job_uuid = firework["spec"]["_tasks"][0]["job"]["uuid"]
                docs = docs_store.query(
                    {"uuid": job_uuid, "output.properties": {"$exists": True}}
                )
                for doc in docs:
                    blob_uuid = doc["output"]["properties"]["blob_uuid"]
                    _properties = property_store.query_one({"blob_uuid": blob_uuid})
                    yield {'properties': _properties, 'cid': cid}
    def process_item(self, item: Dict) -> Dict:
        new_item = dict(**item['properties']['data'])
        new_item['cid'] = item['cid']
        return new_item
    def update_targets(self, items: List[Dict]):
        for item in items:
            target_doc = self.target_store.query_one({'cid': item['cid']})
            if target_doc is not None:
                for field in self.fields:
                    if field in item: 
                        if field not in target_doc:
                            target_doc[field] = {}
                        target_doc[field] = item[field] | target_doc[field]
                self.target_store.update(target_doc, "cid")
            else: 
                self.target_store.update(item, "cid")

#
builder = MolecularProperties(workflow_store, fireworks_store, docs_store, property_store, molecule_property_store)
builder.run(log_level=logging.INFO)
#
