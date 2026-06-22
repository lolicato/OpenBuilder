import json


def run_from_config(config_path: str):
    with open(config_path) as f:
        data = json.load(f)

    for builder_name, config_dict in data.items():
        config_module = __import__(f"src.builders.{builder_name}.config", fromlist=["Config"])
        config = config_module.Config(**config_dict)

        runner_module = __import__(f"src.builders.{builder_name}.runner", fromlist=[""])
        runner_class = next(c for name, c in vars(runner_module).items() if name.endswith("Runner"))
        runner_class().run(config)