import os
import sys
import json
import csv
from pathlib import Path
import subprocess
import multiprocessing
import logging
from datetime import datetime
import shutil
import csv
import glob
import pandas as pd
import operator
import argparse


# explicit paths to chemprop_predict and sklearn_predict
path_to_chemprop_predict = shutil.which("chemprop_predict")
path_to_sklearn_predict = shutil.which("sklearn_predict")
path_to_chemfunc = shutil.which("chemfunc")

level = logging.WARNING
logger = logging.getLogger(__name__)
logger.setLevel(level)

# detailed log
# formatter = logging.Formatter(
#     '%(levelname)s %(asctime)s : (%(filename)s::%(funcName)s::%(lineno)d) : %(message)s')

# basic log
formatter = logging.Formatter('%(levelname)s %(asctime)s : %(message)s')

stream_handler = logging.StreamHandler(sys.stdout)
stream_handler.setFormatter(formatter)

file_handler = logging.FileHandler('eskape-app.log')
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(formatter)

logger.addHandler(stream_handler)
logger.addHandler(file_handler)


def run_random_forest_models(test_path, checkpoint_dir, preds_path):
    logger.info("run_random_forest_models ...")
    try:
        result = subprocess.run(
            f"{path_to_sklearn_predict} --test_path {test_path} --checkpoint_dir {checkpoint_dir} --preds_path {preds_path}", shell=True, check=True, capture_output=True)
        # success
        if result.returncode == 0:
            logger.info("Success")
            return True
        else:
            logger.info("Fail")
            return False
    except subprocess.CalledProcessError as e:
        logger.info("Fail")
        return False


def run_chemprop_models(test_path, checkpoint_dir, preds_path):
    logger.info("run_chemprop_models ...")
    try:
        result = subprocess.run(
            f"{path_to_chemprop_predict} --individual_ensemble_predictions --test_path {test_path} --checkpoint_dir {checkpoint_dir} --preds_path {preds_path}", shell=True, check=True, capture_output=True)
        # success
        if result.returncode == 0:
            logger.info("Success")
            return True
        else:
            logger.info("Fail")
            return False
    except subprocess.CalledProcessError as e:
        logger.info("Fail")
        return False


def run_chemprop_rdkit_models(test_path, checkpoint_dir, preds_path):
    logger.info("run_chemprop_rdkit_models ...")
    try:
        result = subprocess.run(f"{path_to_chemprop_predict} --individual_ensemble_predictions --test_path {test_path} --checkpoint_dir {checkpoint_dir} --preds_path {preds_path} --features_generator rdkit_2d_normalized --no_features_scaling", shell=True, check=True, capture_output=True)

        # success
        if result.returncode == 0:
            logger.info("Success")
            return True
        else:
            logger.info("Fail")
            return False
    except subprocess.CalledProcessError as e:
        logger.info(f"Fail")
        return False


def get_canonical_data(f1, f2, f3):
    # reading two csv files
    data1 = pd.read_csv(f1)
    data2 = pd.read_csv(f2)

    # rename coloumns for data2
    data2.rename(columns={'smiles': 'user_smiles'}, inplace=True)
    data2.rename(columns={'canonical': 'smiles'}, inplace=True)

    # using merge csv on smiles column
    output = pd.merge(data1, data2,  on='smiles', how='inner')

    # output
    output.to_csv(f"{f3}-validated.txt", index=False)


def predict_function(args):
    # info
    logger.info(f"path_to_chemprop_predict: {path_to_chemprop_predict}")
    logger.info(f"path_to_sklearn_predict: {path_to_sklearn_predict}")
    logger.info(f"path_to_chemfunc: {path_to_chemfunc}")

    path, filename = os.path.split(args.input_file)
    # logger.info(f"path: {path}")
    # logger.info(f"filename: {filename}")

    base_directory = args.models_directory
    # logger.info(f"base_directory: {base_directory}")

    try:
        Path(path).mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logger.info(f"Folder {path} is already there: {e}")
    else:
        # logger.info(f"Folder '{path}' was created")
        pass

    # input SMILES to predict on
    test_path = args.input_file

    # path where a CSV file containing the predictions will be saved
    output_directory = args.output_directory

    try:
        Path(output_directory).mkdir(parents=True, exist_ok=False)
    except Exception as e:
        logger.info(f"Folder {output_directory} is already there: {e}")
    else:
        logger.info(f"Folder '{output_directory}' was created")

    preds_path = os.path.join(output_directory, filename)

    # calculate TNN (Tanimoto nearest neighbour) metric
    # TNN script: https://github.com/swansonk14/chemfunc/blob/main/src/chemfunc/nearest_neighbor.py
    reference_path = base_directory+"/canonical_data/training_data_canonical.csv"
    result = subprocess.run(
        f"{path_to_chemfunc} nearest_neighbor --data_path {test_path} --save_path {preds_path}-nearest-neighbor.txt --metric tanimoto --reference_smiles_column smiles --reference_data_path {reference_path} > /dev/null 2>&1", shell=True, check=True)

    if result.returncode == 0:
        logger.info("Success")
    else:
        logger.error(
            f"File {test_path} failed for 'chemfunc nearest_neighbor' funcion")

    # calculate molecular weight and clogP using user submitted smile
    # https://github.com/swansonk14/chemfunc/blob/main/chemfunc/molecular_properties.py
    # get Canonical SMILES using https://github.com/swansonk14/chemfunc/blob/main/chemfunc/canonicalize_smiles.py
    result = subprocess.run(
        f"{path_to_chemfunc} compute_properties --data_path {test_path} --save_path {preds_path}-properties.txt --properties mol_weight clogp > /dev/null 2>&1", shell=True, check=True)

    if result.returncode == 0:
        logger.info("Success")
    else:
        logger.error(
            f"File {test_path} failed for 'chemfunc compute_properties' funcion")

    # get Canonical SMILES using https://github.com/swansonk14/chemfunc/blob/main/chemfunc/canonicalize_smiles.py
    result = subprocess.run(
        f"{path_to_chemfunc} canonicalize_smiles --data_path {test_path} --save_path {preds_path}-canonical.txt --remove_salts --smiles_column smiles --canonical_smiles_column canonical > /dev/null 2>&1", shell=True, check=True)

    if result.returncode == 0:
        logger.info("Success")
        # parse the input canonical and match-up to training canonical smiles
        get_canonical_data(base_directory+"/canonical_data/training_data_canonical.csv",
                           f"{preds_path}-canonical.txt", preds_path)
    else:
        logger.error(
            f"File {test_path} failed for 'chemfunc canonicalize_smiles' funcion")

    # run in parallel
    try:

        results = process_jobs(test_path, base_directory, preds_path)

        if results == 0:
            logger.info("Done.")
            if args.debug == False:
                print("Done.")
        else:
            logger.error(
                f"File {test_path} not processed properly return code {results}")

    except Exception as e:
        print(f"Unexpected {e=}, {type(e)=}")
        logger.error(f"{e}")


def launch_job(test_path, base_directory, preds_path, pathogen):
    logger.info(f"launch_job for {pathogen} ...")
    try:
        p1 = multiprocessing.Process(target=run_chemprop_rdkit_models, args=(
            f"{test_path}", f"{base_directory}/models/all/{pathogen}_rdkit", f"{preds_path}-{pathogen}_rdkit.csv",))
        p2 = multiprocessing.Process(target=run_chemprop_models, args=(
            f"{test_path}", f"{base_directory}/models/all/{pathogen}_chemprop", f"{preds_path}-{pathogen}_chemprop.csv",))
        p3 = multiprocessing.Process(target=run_random_forest_models, args=(
            f"{test_path}", f"{base_directory}/models/all/{pathogen}_rf", f"{preds_path}-{pathogen}_rf.csv",))
        p1.run()
        p2.run()
        p3.run()
        logger.info("Running threads ...")
    except Exception as e:
        print(e)
        logger.error(f"{e}")


def worker(test_path, base_directory, preds_path):
    launch_job(test_path, base_directory, preds_path, "AB")
    launch_job(test_path, base_directory, preds_path, "BW")
    launch_job(test_path, base_directory, preds_path, "DKO")
    launch_job(test_path, base_directory, preds_path, "EF")
    launch_job(test_path, base_directory, preds_path, "KP")
    launch_job(test_path, base_directory, preds_path, "PA")
    launch_job(test_path, base_directory, preds_path, "SA")


def process_jobs(test_path, base_directory, preds_path):
    try:
        main_thread = multiprocessing.Process(target=worker, args=(
            f"{test_path}", f"{base_directory}", f"{preds_path}",))
        main_thread.run()
        main_thread.join()
        return main_thread.exitcode
    except Exception as e:
        if main_thread.exitcode is None:
            parse_results(f"{preds_path}")
            logger.info("analysis complete.")
            return 0
        else:
            # exit unsuccessfully
            logger.error("analysis failed.")
            print(f"exit unsuccessfully, {e}")
            return 1


def get_properties(o_f_path, accession, smile):
    with open(os.path.join(o_f_path, f"{accession}-properties.txt"), "r") as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in reader:
            if row[0] != "smiles":
                if row[0] == smile:
                    return {"mol_weight": str_to_float(row[1]), "clogp": str_to_float(row[2])}


def get_nearest_neighbor(o_f_path, accession, smile):
    with open(os.path.join(o_f_path, f"{accession}-nearest-neighbor.txt"), "r") as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in reader:
            if row[0] != "smiles":
                if row[0] == smile:
                    return {"tanimoto_nearest_neighbor": row[1], "tanimoto_nearest_neighbor_similarity": str_to_float(row[2])}


def str_to_float(input_string):
    converted = 0.0
    try:
        converted = float(input_string)
    except ValueError:
        converted = 0.0
    return converted


def get_validated(o_f_path, accession, smile):
    with open(os.path.join(o_f_path, f"{accession}-validated.txt"), "r") as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in reader:
            if row[0] != "smiles":
                if row[8] == smile:
                    return {"canonical_smiles": row[0],
                            "EF": str_to_float(row[1]),
                            "SA": str_to_float(row[2]),
                            "KP": str_to_float(row[3]),
                            "AB": str_to_float(row[4]),
                            "PA": str_to_float(row[5]),
                            "BW": str_to_float(row[6]),
                            "DKO": str_to_float(row[7]),
                            "user_smiles": row[8],
                            }


def get_sum_ppfs(smile, j):
    results = {}
    for i in j:
        results[i] = {}
        for model_type in j[i][smile]:
            if model_type in ["rf", "chemprop", "rdkit"]:
                top_score = dict(sorted(j[i][smile][model_type].items(
                ), key=operator.itemgetter(1), reverse=True)[:7])
                l = list(top_score.values())
                [str_to_float(x) for x in l]
                results[i][f"sum_{model_type}"] = sum(l)
                try:
                    results[i][f"ppf_{model_type}"] = l[0] / l[1]
                except ZeroDivisionError as e:
                    results[i][f"ppf_{model_type}"] = 0.0

    return results


def filter_results_tabular(filename=False, cplog=False, molweight=False, sum=False, ppf=False, modeltyperdkit=False, modeltypechemprop=False, modeltyperf=False, validated=False, tnn=False):
    # read in tsv and drop columns not needed
    df = pd.read_table(filename)

    if cplog == 'true':
        df.drop(columns=["clogp"], axis=1, inplace=True)
    if tnn == 'true':
        df.drop(columns=["tanimoto_nearest_neighbor_similarity", "tanimoto_nearest_neighbor"],
                axis=1, inplace=True)
    if molweight == 'true':
        df.drop(columns=["molecular_weight"], axis=1, inplace=True)
    if sum == 'true':
        if "sum_rdkit" in df.columns:
            df.drop(columns=["sum_rdkit"], axis=1, inplace=True)
        if "sum_chemprop" in df.columns:
            df.drop(columns=["sum_chemprop"], axis=1, inplace=True)
        if "sum_rf" in df.columns:
            df.drop(columns=["sum_rf"], axis=1, inplace=True)

    if ppf == 'true':
        if "ppf_rdkit" in df.columns:
            df.drop(columns=["ppf_rdkit"], axis=1, inplace=True)
        if "ppf_chemprop" in df.columns:
            df.drop(columns=["ppf_chemprop"], axis=1, inplace=True)
        if "ppf_rf" in df.columns:
            df.drop(columns=["ppf_rf"], axis=1, inplace=True)

    if modeltyperdkit == 'true':
        # delete all rdkit values for each pathogen
        df.drop(columns=["ef_rdkit", "sa_rdkit", "kp_rdkit",
                         "ab_rdkit", "pa_rdkit", "bw_rdkit", "dko_rdkit"], axis=1, inplace=True)
        if "sum_rdkit" in df.columns:
            df.drop(columns=["sum_rdkit"], axis=1, inplace=True)
        if "ppf_rdkit" in df.columns:
            df.drop(columns=["ppf_rdkit"], axis=1, inplace=True)

    if modeltypechemprop == 'true':
        # delete all chemprop values for each pathogen
        df.drop(columns=["ef_chemprop", "sa_chemprop", "kp_chemprop",
                "ab_chemprop", "pa_chemprop", "bw_chemprop", "dko_chemprop"], axis=1, inplace=True)
        if "sum_chemprop" in df.columns:
            df.drop(columns=["sum_chemprop"], axis=1, inplace=True)
        if "ppf_chemprop" in df.columns:
            df.drop(columns=["ppf_chemprop"], axis=1, inplace=True)

    if modeltyperf == 'true':
        # delete all rf values for each pathogen
        df.drop(columns=["ef_rf", "sa_rf", "kp_rf",
                "ab_rf", "pa_rf", "bw_rf", "dko_rf"], axis=1, inplace=True)
        if "sum_rf" in df.columns:
            df.drop(columns=["sum_rf"], axis=1, inplace=True)
        if "ppf_rf" in df.columns:
            df.drop(columns=["ppf_rf"], axis=1, inplace=True)

    if validated == 'true':
        # delete all rf values for each pathogen
        df.drop(columns=["ef_validated", "sa_validated", "kp_validated",
                "ab_validated", "pa_validated", "bw_validated", "dko_validated"], axis=1, inplace=True)

    # write tabular
    o_f_path, f_name = os.path.split(filename)
    # basename
    fn = f_name.split("-final.tsv")[0]
    # logger.info("save modified tab %s", os.path.join(
    #     o_f_path, f"{fn}-final-mod.tsv"))

    df.to_csv(os.path.join(o_f_path, f"{fn}-final-mod.tsv"),
              sep='\t', index=False, header=True)


def generate_results_tabular(o_f_path, accession):
    j = ""
    try:
        with open(os.path.join(o_f_path, f"{accession}-final.json"), 'r') as jfile:
            j = json.load(jfile)
    except Exception as e:
        print(e)
        logger.error("failed to generate results, %s", e)
        exit()

    with open(os.path.join(o_f_path, f"{accession}-final.tsv"), "w") as tabout:
        out = {}
        writer = csv.writer(tabout, delimiter='\t', dialect='excel')
        writer.writerow(["smiles", "molecular_weight", "clogp",
                        "sum_rdkit", "sum_chemprop", "sum_rf",
                         "ppf_rdkit", "ppf_chemprop", "ppf_rf",
                         "ef_rdkit", "ef_chemprop", "ef_rf",
                         "sa_rdkit", "sa_chemprop", "sa_rf",
                         "kp_rdkit", "kp_chemprop", "kp_rf",
                         "ab_rdkit", "ab_chemprop", "ab_rf",
                         "pa_rdkit", "pa_chemprop", "pa_rf",
                         "bw_rdkit",  "bw_chemprop", "bw_rf",
                         "dko_rdkit", "dko_chemprop", "dko_rf",
                         "ef_validated", "sa_validated", "kp_validated",
                         "ab_validated", "pa_validated", "bw_validated",
                         "dko_validated", "tanimoto_nearest_neighbor", "tanimoto_nearest_neighbor_similarity"
                         ])
        for i in j:
            for smile in j[i]:
                out[smile] = [smile, j[i][smile]["molweight"], j[i][smile]["clogp"],
                              j[i][smile]["sum_rdkit"], j[i][smile]["sum_chemprop"], j[i][smile]["sum_rf"],
                              j[i][smile]["ppf_rdkit"], j[i][smile]["ppf_chemprop"], j[i][smile]["ppf_rf"],
                              j[i][smile]["rdkit"]["EF"], j[i][smile]["chemprop"]["EF"], j[i][smile]["rf"]["EF"],
                              j[i][smile]["rdkit"]["SA"], j[i][smile]["chemprop"]["SA"], j[i][smile]["rf"]["SA"],
                              j[i][smile]["rdkit"]["KP"], j[i][smile]["chemprop"]["KP"], j[i][smile]["rf"]["KP"],
                              j[i][smile]["rdkit"]["AB"], j[i][smile]["chemprop"]["AB"], j[i][smile]["rf"]["AB"],
                              j[i][smile]["rdkit"]["PA"], j[i][smile]["chemprop"]["PA"], j[i][smile]["rf"]["PA"],
                              j[i][smile]["rdkit"]["BW"], j[i][smile]["chemprop"]["BW"], j[i][smile]["rf"]["BW"],
                              j[i][smile]["rdkit"]["DKO"], j[i][smile]["chemprop"]["DKO"], j[i][smile]["rf"]["DKO"],
                              j[i][smile]["validated"]["EF"], j[i][smile]["validated"]["SA"], j[i][smile]["validated"]["KP"],
                              j[i][smile]["validated"]["AB"], j[i][smile]["validated"]["PA"], j[i][smile]["validated"]["BW"],
                              j[i][smile]["validated"]["DKO"], j[i][smile]["tanimoto_nearest_neighbor"], j[i][smile]["tanimoto_nearest_neighbor_similarity"]
                              ]
        for key, value in out.items():
            writer.writerow(value)


def parse_results(path):
    logger.info(f"parse_results at {path} ...")
    """
    Read the csvs and create a json file with all the results

    Format:
    'accession': { 
            'smile_1': { 
                rdkit: {EF: val_1, AB: val_1, ..} ,
                chemprop: {EF: val_1, AB: val_1, ..} ,
                rf: {EF: val_1, AB: val_1, ..} ,
                validated: {EF: val_1, AB: val_1, ..},
                molweight: "val_1",
                clogp: "val_1",
                tanimoto_nearest_neighbor: "val_1",
                tanimoto_nearest_neighbor_similarity": "val_1",
                sumofpredictionscores: "val_1",
                ppf: "val_1"

            },
            'smile_2': { 
                rdkit: {EF: val_2, AB: val_2, ..} ,
                chemprop: {EF: val_2, AB: val_2, ..} ,
                rf: {EF: val_2, AB: val_2, ..} ,
                validated: {EF: val_2, AB: val_2, ..} ,
                molweight: "val_2",
                clogp: "val_2",
                tanimoto_nearest_neighbor: "val_2",
                tanimoto_nearest_neighbor_similarity": "val_2",
                sumofpredictionscores: "val_2",
                ppf: "val_2"
            },
            'smile_N': { 
                rdkit: {EF: val_N, AB: val_N, ..} ,
                chemprop: {EF: val_N, AB: val_N, ..} ,
                rf: {EF: val_N, AB: val_N, ..} ,
                validated: {EF: val_N, AB: val_N, ..},
                molweight: "val_N",
                clogp: "val_N",
                tanimoto_nearest_neighbor: "val_N",
                tanimoto_nearest_neighbor_similarity": "val_N",
                sumofpredictionscores: "val_N",
                ppf: "val_N"
            },
            ...
    }
    """
    o_f_path, o_f_name = os.path.split(path)
    accession = o_f_name
    results = {}
    results[accession] = {}
    files = glob.glob(os.path.join(o_f_path, "*"))
    for f in files:
        if os.path.isfile(f) and os.path.splitext(os.path.basename(f))[1][1:].strip() in ["csv"]:
            with open(f, 'r') as csvfile:
                reader = csv.reader(csvfile, delimiter=',', quotechar='|')
                filename = os.path.splitext(os.path.basename(f))[0]
                _filename = filename.split("-")[-1].split("_")
                pathogen = _filename[0]
                model_type = _filename[1]
                for row in reader:
                    if "smiles" not in row[0]:
                        if len(results[accession]) == 0:
                            results[accession] = {
                                row[0]: {
                                    "rf": {"DKO": "", "SA": "", "PA": "", "EF": "", "BW": "", "KP": "", "AB": ""},
                                    "chemprop": {"DKO": "", "SA": "", "PA": "", "EF": "", "BW": "", "KP": "", "AB": ""},
                                    "rdkit": {"DKO": "", "SA": "", "PA": "", "EF": "", "BW": "", "KP": "", "AB": ""},
                                    "validated": {"DKO": "", "SA": "", "PA": "", "EF": "", "BW": "", "KP": "", "AB": ""},
                                    "molweight": "",
                                    "clogp": "",
                                    "tanimoto_nearest_neighbor": "",
                                    "tanimoto_nearest_neighbor_similarity": ""
                                }
                            }
                        elif row[0] not in results[accession].keys():
                            results[accession].update({row[0]: {
                                "rf": {"DKO": "", "SA": "", "PA": "", "EF": "", "BW": "", "KP": "", "AB": ""},
                                "chemprop": {"DKO": "", "SA": "", "PA": "", "EF": "", "BW": "", "KP": "", "AB": ""},
                                "rdkit": {"DKO": "", "SA": "", "PA": "", "EF": "", "BW": "", "KP": "", "AB": ""},
                                "validated": {"DKO": "", "SA": "", "PA": "", "EF": "", "BW": "", "KP": "", "AB": ""},
                                "molweight": "",
                                "clogp": "",
                                "tanimoto_nearest_neighbor": "",
                                "tanimoto_nearest_neighbor_similarity": ""
                            }})

                        results[accession][row[0]
                                           ][model_type][pathogen] = str_to_float(row[1])

    # write json file with results per session
    with open(os.path.join(o_f_path, f"{accession}.json"), "w") as outfile:
        json.dump(results, indent=2, fp=outfile)

    # read the results json file
    try:
        with open(os.path.join(o_f_path, f"{accession}.json"), 'r') as jsonfile:
            j = json.load(jsonfile)
    except Exception as e:
        print(e)
        exit()

    # add validated, molweight, clogp, sumofpredictionscores, and ppf to the results
    properties = {}
    nearest_neighbor = {}
    validated = {}
    sums_ppfs = {}

    for i in j:
        #  for each smile, get properties
        for k in j[i].keys():
            properties[k] = get_properties(o_f_path, accession, k)
            nearest_neighbor[k] = get_nearest_neighbor(o_f_path, accession, k)
            validated[k] = get_validated(o_f_path, accession, k)
            sums_ppfs[k] = get_sum_ppfs(k, j)

    # update the json file
    for index in j:
        for m in j[index].keys():
            # add properties
            if properties[m]:
                j[index][m]["molweight"] = properties[m]["mol_weight"]
                j[index][m]["clogp"] = properties[m]["clogp"]
            # add nearest_neighbor
            if nearest_neighbor[m]:
                j[index][m]["tanimoto_nearest_neighbor"] = nearest_neighbor[m]["tanimoto_nearest_neighbor"]
                j[index][m]["tanimoto_nearest_neighbor_similarity"] = nearest_neighbor[m]["tanimoto_nearest_neighbor_similarity"]
            # add validated
            if validated[m]:
                for x in validated[m]:
                    j[index][m]["validated"][x] = validated[m][x]
            # add sums
            if sums_ppfs[m]:
                j[index][m]["sum_rdkit"] = sums_ppfs[m][index]["sum_rdkit"]
                j[index][m]["sum_chemprop"] = sums_ppfs[m][index]["sum_chemprop"]
                j[index][m]["sum_rf"] = sums_ppfs[m][index]["sum_rf"]
            # add ppf
            if sums_ppfs[m]:
                j[index][m]["ppf_rdkit"] = sums_ppfs[m][index]["ppf_rdkit"]
                j[index][m]["ppf_chemprop"] = sums_ppfs[m][index]["ppf_chemprop"]
                j[index][m]["ppf_rf"] = sums_ppfs[m][index]["ppf_rf"]

    # write json file with results per session
    with open(os.path.join(o_f_path, f"{accession}-final.json"), "w") as outfile1:
        json.dump(j, indent=2, fp=outfile1)

    # write tabular file with results
    generate_results_tabular(o_f_path, accession)


def main():
    parser = argparse.ArgumentParser(
        prog="eskape_model", description="standalone")
    parser.add_argument('-i', '--input_file', dest="input_file",
                        help="input csv file with SMILES, with header 'smiles'")
    parser.add_argument('-m', '--models_directory', dest="models_directory",
                        help="path to trained models")
    parser.add_argument('-o', '--output_directory', dest="output_directory",
                        help="path to save results")
    parser.add_argument('--debug', dest="debug",
                        action="store_true", help="debug mode")

    if len(sys.argv) == 1:
        sys.stderr.write("No arguments provided, printing help menu ...\n")
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    if args.debug:
        logger.setLevel(10)

    predict_function(args)


if __name__ == '__main__':
    main()
