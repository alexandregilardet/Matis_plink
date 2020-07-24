#!/usr/bin/env python3

import sys
import pandas as pd
from copy import deepcopy

# READ INPUT


def reader(file):
    """reads and appends lines of a file to a list"""
    data_lines = []
    with open(file, "r") as f:
        for line in f.readlines():
            data_lines.append(line)
    return data_lines[:]


def from_ped(data_lines_ped):
    """parses the list of lines from .ped file into a list of dict"""
    parsed_ped = []
    parsed_line = {}
    for line in data_lines_ped:
        line = line.strip().split("\t")
        parsed_line["family"] = line[0]
        parsed_line["individual_ID"] = line[1]
        parsed_line["parent_A"] = line[2]
        parsed_line["parent_B"] = line[3]
        parsed_line["sex"] = line[4]
        parsed_line["phenotype"] = line[5]
        parsed_line["loci"] = line[6:]
        parsed_ped.append(parsed_line.copy())
        # copy to avoid reference to dict which iterates
    return parsed_ped[:]


def from_map(data_lines_map):
    """parses the list of lines from .map file into a list of dict"""
    parsed_map = []
    parsed_line = {}
    for line in data_lines_map:
        line = line.strip().split("\t")
        parsed_line["chromosome"] = line[0]
        parsed_line["variant_ID"] = line[1]
        parsed_line["position"] = line[2]
        parsed_line["coordinate"] = line[3]
        parsed_map.append(parsed_line.copy())
    return parsed_map[:]


def from_ped_map(parsed_ped, parsed_map):
    """fuses alleles (from .ped) and other locus info (from .map) into a common dictionnary attributed
    to each individual and for each marker"""
    for ind in parsed_ped:
        marker_loci = []
        for locus, marker in zip(ind["loci"], parsed_map):
            marker["alleles"] = locus
            marker_loci.append(marker)
        marker_loci_copy = deepcopy(marker_loci)  # not to affect the original dict
        ind["loci"] = marker_loci_copy
    return parsed_ped[:]  # which is now a fusion with parsed_map


def check_files(parsed_ped, parsed_map):
    """checks if .ped and .map files are likely to correspond to same dataset"""
    try:
        nb_loci = len(parsed_ped[1]["loci"])
        nb_marker = len(parsed_map)
        assert nb_loci == nb_marker
    except AssertionError:
        sys.exit(
            f"""Something wrong with number of markers, {nb_loci} loci in .ped file and 
        {nb_marker} variants in .map file. Files might not be from same data."""
        )
    return None


# FILTERING


def check_target(target):
    """asserts that target is not a string and does not contain duplicates"""
    try:
        if isinstance(target, str):
            sys.exit("You provided a string as target instead of a list of strings")
        assert len(set(target)) - len(target) == 0
    except AssertionError:
        sys.exit("There are duplicates in your query")
    return None


def get_individuals(target, parsed_ped_map):
    """extracts a subset of individuals and returns an object similar to parsed_ped_map"""
    try:
        check_target(target)
        found = []  # to check if everything was found
        new_ind = []
        for elmt in target:
            for ind in parsed_ped_map:
                if ind["individual_ID"] == elmt:
                    new_ind.append(ind)
                    found.append(elmt)
        diff = list(set(target) - set(found))
        if diff != []:
            print(
                f"The individual(s) {diff} from target could not be found in available data"
            )
    except TypeError:
        sys.exit("Target should be a list of strings")
    return new_ind[:]


def get_markers(target, parsed_ped_map):
    """extracts a subset of markers and returns an object similar to parsed_ped_map"""
    try:
        check_target(target)
        copy_dict = deepcopy(parsed_ped_map)  # avoid modifying original
        found = []
        for ind in copy_dict:
            new_markers = []
            for marker in ind["loci"]:
                for elmt in target:
                    if marker["variant_ID"] == elmt:
                        new_markers.append(marker)
                        found.append(elmt)
            ind["loci"] = new_markers
        diff = list(set(target) - set(found))
        if diff != []:
            print(
                f"The marker(s) {diff} from target could not be found in available data"
            )
    except TypeError:
        sys.exit("Target should be a list of strings")
    return copy_dict[:]


# EXPORT


def from_fused_to_ped_map(parsed_ped_map):
    """parses the parsed fused object back to both map and ped parsed objects"""
    pedkey = ["alleles"]
    mapkeys = ["chromosome", "coordinate", "position", "variant_ID"]
    new_map = []
    copy_dict2 = deepcopy(parsed_ped_map)
    filterByKey = lambda keys: {x: marker[x] for x in keys}
    for marker in copy_dict2[0]["loci"]:  # new map
        map_marker = filterByKey(mapkeys)
        new_map.append(map_marker)
    for ind in copy_dict2:  # new ped
        new_ped = []
        for marker in ind["loci"]:
            ped = filterByKey(pedkey)
            new_ped.append(ped)
        gt = []
        for alleles in new_ped:
            allele = list(alleles.values())
            gt += allele
        ind["loci"] = gt
    return new_map[:], copy_dict2[:]


def to_map(new_map):
    """parses the map object back to map lines that can be further saved in a .map file"""
    map_lines_to_save = []
    for marker in new_map:
        line = list(marker.values())
        line[1], line[3] = line[3], line[1]  # swap order of coordinate and variant_ID
        line = "\t".join(line)
        line = line + "\n"
        map_lines_to_save.append(line)
    return map_lines_to_save[:]


def to_ped(copy_dict2):
    """parses the ped object back to ped lines that can be further saved in a .ped file"""
    ped_lines_to_save = []
    for ind in copy_dict2:
        line = list(ind.values())
        line += line[6]
        line.remove(line[6])  # to end with alleles
        line = "\t".join(line)
        line = line + "\n"
        ped_lines_to_save.append(line)
    return ped_lines_to_save[:]


def writer(content, name):
    """writes a file"""
    with open(name, "w") as f:
        f.writelines(content)
    return None


def export_to_ped_map(parsed_ped_map, name):
    """main function to parse back a parsed fused object and export it as .ped and .map files"""
    new_map, copy_dict2 = from_fused_to_ped_map(parsed_ped_map)
    map_lines_to_save = to_map(new_map)
    ped_lines_to_save = to_ped(copy_dict2)
    name_map = name + ".map"
    writer(map_lines_to_save, name_map)
    name_ped = name + ".ped"
    writer(ped_lines_to_save, name_ped)
    return None


def from_fused_to_df(parsed_ped_map):
    """returns a data frame as individual x marker = genotype from the parsed fused object"""
    tbl = []
    for ind in parsed_ped_map:
        line = {}
        for marker in ind["loci"]:
            line["ind"] = ind["individual_ID"]
            line["marker"] = marker["variant_ID"]
            line["gt"] = marker["alleles"]
            line_copy = deepcopy(line)
            tbl.append(line_copy)
    df = pd.DataFrame(tbl)
    return df[:]


def reshape_df(df):  # order of markers changed compared to original .map file
    """reshapes a data frame from long to wide format"""
    df2 = df.pivot_table(  # long to wide format
        index="ind",
        columns="marker",
        values="gt",  # str so need to adjust aggfunc
        aggfunc=lambda x: " ".join(
            x
        ),  # if several values for same index, will concatenate with space
    )
    return df2[:]


def writer_to_csv(df2, name):
    """writes data frame to csv file"""
    file = name + ".csv"
    df2.to_csv(file, sep=";")
    return None


def export_to_csv(parsed_ped_map, name):
    """main function to export parsed fused object to csv file"""
    df = from_fused_to_df(parsed_ped_map)
    df2 = reshape_df(df)
    writer_to_csv(df2, name)
    return None


def gt_to_alleles(gt):
    """parses plink genotypes into both alleles"""
    allele_1 = gt.split(" ")[0]
    allele_2 = gt.split(" ")[1]
    return allele_1[:], allele_2[:]


def allele_to_genepop(allele):
    """translates plink allele into genepop format"""
    if allele == "A":
        allele_gen = "001"
    elif allele == "B":
        allele_gen = "002"
    elif allele == "0":
        allele_gen = "000"
    return allele_gen[:]


def allele_to_structure(allele):
    """translates plink allele into structure format"""
    if allele == "A":
        allele_str = "1"
    elif allele == "B":
        allele_str = "2"
    elif allele == "0":
        allele_str = "9"
    return allele_str[:]


def gt_parser(df, format):
    """parses and translates plink genotypes into either genepop or structure format"""
    new_gt = []
    for gt in df["gt"]:
        allele_1, allele_2 = gt_to_alleles(gt)
        if format == "genepop":
            new_1 = allele_to_genepop(allele_1)
            new_2 = allele_to_genepop(allele_2)
        elif format == "structure":
            new_1 = allele_to_structure(allele_1)
            new_2 = allele_to_structure(allele_2)
        new_gt.append(new_1 + new_2)
    df["gt"] = new_gt
    return df[:]


def parser_to_genepop(df2):
    """parses wide format data frame to genepop format"""
    df2.insert(0, "", ",")  # add column with comma
    header = []
    for col in df2.columns:
        col = col + "\n"
        header.append(col)
    header.pop(0)  # remove first element '' from the comma column
    return df2[:], header[:]


def writer_to_genepop(df3, header, name, title):
    """writes a genepop file"""
    with open(name, "w") as f:
        f.write('Title line: "' + title + '"\n')
        f.writelines(header)
        f.write("Pop\n")
        df3.to_csv(f, header=False, sep=" ")
    return None


def export_to_genepop(parsed_ped_map, name):
    """main function to export parsed fused object to genepop format"""
    df = from_fused_to_df(parsed_ped_map)
    df_allele = gt_parser(df, "genepop")
    df2 = reshape_df(df_allele)
    df3, header = parser_to_genepop(df2)
    title = input("Enter the title line for genepop file:")
    file = name + ".txt"
    writer_to_genepop(df3, header, file, title)
    return None


def reshape_to_structure(df):
    """reshapes data frame on two columns per marker for structure format"""
    header = ""  # create header line with marker once instead of marker.1 marker.2 like in df
    for col in df.columns[0:]:
        df["{}.1".format(col)], df["{}.2".format(col)] = zip(
            *df[col].apply(lambda x: list(x))
        )
        del df[col]  # remove marker on only 1 column
        header += col
        header += " "
    header.rstrip()  # remove last blank
    return df[:], header[:]


def writer_to_structure(df3, header, file):
    """writes data frame to structure file"""
    with open(file, "w") as f:
        f.write(header + "\n")
        df3.to_csv(f, header=False, sep=" ")
    return None


def export_to_structure(parsed_ped_map, name):
    """main function to export parsed fused object to structure format"""
    df = from_fused_to_df(parsed_ped_map)
    df_allele = gt_parser(df, "structure")
    df2 = reshape_df(df_allele)
    df3, header = reshape_to_structure(df2)
    file = name + ".str"
    writer_to_structure(df3, header, file)
    return None


def export(obj, format, name):
    """general function to call the right export function according to the format specified"""
    try:
        assert isinstance(name, str)
        if format == "ped_map":
            export_to_ped_map(obj, name)
        elif format == "csv":
            export_to_csv(obj, name)
        elif format == "genepop":
            export_to_genepop(obj, name)
        elif format == "structure":
            export_to_structure(obj, name)
        else:
            sys.exit(f"Called an unsupported export format: {format}")
    except AssertionError:
        sys.exit("File name should be a string")
    except TypeError:
        sys.exit("Tried to export unsupported object")
    return None


class plinkPy:
    """Class whose instance contains related .ped and .map files exported from Plink"""

    def __init__(self, ped, map, parsed=None):  # self instance #file attribute
        self.ped = ped
        self.map = map
        self.parsed = parsed
        # parsed will be the fused parsed_ped_map object
        return None

    def __str__(self):  # for end-user #str() or print()
        return f".ped file {self.ped} - .map file {self.map}"

    def __repr__(self):  # for other developer #repr()
        return f".ped file {self.ped} - .map file {self.map}"

    def read(self):
        """main function to read and parse .ped and .map files"""
        try:
            ped_lines = reader(self.ped)
            map_lines = reader(self.map)
            parsed_ped = from_ped(ped_lines)
            parsed_map = from_map(map_lines)
            check_files(parsed_ped, parsed_map)
            parsed_ped_map = from_ped_map(parsed_ped, parsed_map)
        except IndexError:
            sys.exit(
                "Check if .ped and .map files have been inverted, e.g. print your plinkPy object"
            )
        return parsed_ped_map[:]

    @property  # method accessed as attribute
    def individuals(self):
        """prints all individuals and their number"""
        count = 0
        ped_lines = reader(self.ped)
        parsed_ped = from_ped(ped_lines)
        for line in parsed_ped:
            count += 1
            print(line["individual_ID"])
        print(f"There are {count} individuals")
        return None

    @property
    def markers(self):
        """prints all markers and their number"""
        count = 0
        map_lines = reader(self.map)
        parsed_map = from_map(map_lines)
        for line in parsed_map:
            count += 1
            print(line["variant_ID"])
        print(f"There are {count} markers")
        return None

