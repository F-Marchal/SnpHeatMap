import main
import os
import json
import hashlib

def calculate_checksum(file_path, algorithm="sha256", buffer = 4096):
    hash_function = hashlib.new(algorithm)

    with open(file_path, "rb") as flux:
        for chunk in iter(lambda: flux.read(buffer), b""):
            hash_function.update(chunk)

    return hash_function.hexdigest()


def equal_checksums(path1, path2, algorithm="sha256"):
    checksum1 = calculate_checksum(path1, algorithm)
    checksum2 = calculate_checksum(path2, algorithm)
    return checksum1 == checksum2


def tsv_test(name: str, command: str):
    main.main_using_getopts(command + f" -j {name}")
    return compare_results(f"testOutput/{name}/", f"expectedResults/{name}/")


def replace_by_parent_path_wen_needed(path: str, parent: str, sep: str = "/", use_path_as_prefix: bool = True) -> str:
    parent = parent.replace("\\", "/")
    if path[:2] == "./":
        start = -2 if use_path_as_prefix else 0
        return "/".join(parent.split("/")[start:-1]) + sep + path[2:]
    elif path[:3] == "../":
        start = -3 if use_path_as_prefix else 0
        return "/".join(parent.split("/")[start:-2]) + sep + path[2:]
    return path


def compare_results(result: str, expected: str):
    if result[-1] not in ("/", "\\"):
        result += "/"

    if expected[-1] not in ("/", "\\"):
        expected += "/"

    if not os.path.exists(expected + ".test.json"):
        raise FileNotFoundError(f"Can not start a test in {expected}. Can not find .test.json")

    with open(expected + ".test.json", "r", encoding="UTF-8") as json_flux:
        expected_results = json.load(json_flux)

        for path, value in expected_results.copy().items():
            n_path = replace_by_parent_path_wen_needed(path, expected, "_", use_path_as_prefix=True)
            value = replace_by_parent_path_wen_needed(value, expected, "/", use_path_as_prefix=False)

            del expected_results[path]
            expected_results[n_path] = value

    result_list = os.listdir(result)

    logs = {"Unexpected_files": set(),
            "Missing_files": set(result_list).difference(result_list),
            "Nonequivalent": set()}

    for files in result_list:
        if files not in expected_results:
            logs["Unexpected_files"].add(files)
            continue

        if not equal_checksums(expected_results[files], result + files):
            logs["Nonequivalent"].add(files)

    return logs


def main_test():
    test_list = [
        (tsv_test, "Tsv-heatmap_global", "Contig_name BiAllelic_SNP -m 20 -o testOutput -e -1 -wt"),
        (tsv_test, "Tsv-heatmap_global_simplified", "Contig_name BiAllelic_SNP -m 20 -o testOutput -e -1 -wt -i"),
    ]

    final_logs = {}
    for function_, *args in test_list:
        logs = function_(*args)
        name = args[0]
        command = args[1]
        print("Testing " + name + " :", end="")

        for keys, values in logs.copy().items():
            if not values:
                del logs[keys]
            else:
                logs[keys] = list(values)

        if logs:
            logs["Command"] = command
            final_logs[name] = logs
            print("   Failed")

        else:
            print("   Success")

    with open("last_test.json", "w") as json_file:
        json.dump(final_logs, json_file, indent=4)
        print(final_logs)

    return final_logs


main_test()
