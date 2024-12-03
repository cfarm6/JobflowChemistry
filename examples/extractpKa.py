#
from pymongo import MongoClient

#
client = MongoClient(host="localhost", port=27017)
workflows = client["fireworks"]["workflows"]
fireworks = client["fireworks"]["fireworks"]
launches = client["fireworks"]["launches"]


def get_pKa(cid):
    workflow = workflows.find_one({"name": str(cid)})
    nodes = workflow["nodes"]
    firework = fireworks.find_one(
        {"name": "QupKake pKa Prediction", "fw_id": {"$in": nodes}}
    )
    launch_id = firework["launches"][-1]
    launch = launches.find_one({"launch_id": launch_id})
    pKa = launch["action"]["stored_data"]["global"]["pKa"]
    return pKa


df["pKa"] = df["cid"].apply(get_pKa)
df.to_csv("molecules.csv", index=False)
