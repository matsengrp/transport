import json
import os
CONFIG = json.load(open("config.json", "r"))
for key, value in CONFIG.items():
    is_output = "OUTPUT" in key
    does_not_exist = not os.path.exists(str(value))
    if is_output and does_not_exist:
        os.mkdir(value)
