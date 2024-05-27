"""APNS Compressor component

Contributor: pxlxingliang
Date: 2023/12/27
"""

import os
import shutil

def unpack(filepath: str, output_path: str, filetype: str = None, get_support_filetype : bool = False) -> list|str:
    """Unpack file to output_path

    Args:
        filepath (str): file path to decompress
        output_path (str): output path
        filetype (str, optional): specify type of file. Defaults to None.
        get_support_filetype (bool, optional): get supported file types of this function. Defaults to False.

    Raises:
        Exception: if file not exists

    Returns:
        str: output path
        list: supported file types
    """
    # if file type is not supportted, just copy the file to output_path
    if get_support_filetype:
        return ["zip", "tar", "gz", "bz2", "tgz"]
    import zipfile
    import tarfile
    import gzip
    import bz2
    if not os.path.isfile(filepath):
        raise Exception("File %s not exists!" % filepath)
    if not os.path.isdir(output_path):
        os.makedirs(output_path, exist_ok=True)
    if filetype == None:
        if filepath.endswith(".tar.gz"):
            filetype = "tgz"
        else:
            filetype = os.path.splitext(filepath)[1][1:]
    if filetype == "zip":
        with zipfile.ZipFile(filepath, 'r') as zip_ref:
            zip_ref.extractall(output_path)
    elif filetype == "tar":
        with tarfile.open(filepath, "r") as tar_ref:
            tar_ref.extractall(output_path)
    elif filetype == "gz":
        with gzip.open(filepath, 'rb') as gz_ref:
            with open(output_path, 'wb') as out_ref:
                out_ref.write(gz_ref.read())
    elif filetype == "bz2":
        with bz2.open(filepath, 'rb') as bz2_ref:
            with open(output_path, 'wb') as out_ref:
                out_ref.write(bz2_ref.read())
    elif filetype == "tgz":
        with tarfile.open(filepath, "r:gz") as tar_ref:
            tar_ref.extractall(output_path)
    else:
        print("File type \'%s\' is not supported!" % filetype)
        filename = os.path.split(filepath)[-1]
        if os.path.isfile(os.path.join(output_path, filename)):
            os.remove(os.path.join(output_path, filename))
        shutil.copyfile(filepath, os.path.join(output_path, filename))
    print("Unpack %s to %s" % (filepath, output_path))
    return output_path

def pack(packfile_list: list, packfile_name: str, pack_type: str ="zip"):
    """Pack files to packfile_name

    Args:
        packfile_list (list): list of files to pack
        packfile_name (str): name of packed file
        pack_type (str, optional): type of packed file. Defaults to "zip".

    Returns:
        str: name of packed file
    """
    if pack_type == "zip":
        import zipfile
        with zipfile.ZipFile(packfile_name, 'w') as zip_ref:
            for ifile in packfile_list:
                for root, __, ifile in os.walk(ifile):
                    for i in ifile:
                        zip_ref.write(os.path.join(root, i))
    elif pack_type == "tar":
        import tarfile
        with tarfile.open(packfile_name, "w") as tar_ref:
            for ifile in packfile_list:
                tar_ref.add(ifile)
    elif pack_type == "gz":
        import gzip
        with open(packfile_name, 'wb') as out_ref:
            with gzip.open(packfile_name, 'wb') as gz_ref:
                gz_ref.write(out_ref.read())
    elif pack_type == "bz2":
        import bz2
        with open(packfile_name, 'wb') as out_ref:
            with bz2.open(packfile_name, 'wb') as bz2_ref:
                bz2_ref.write(out_ref.read())
    elif pack_type == "tgz":
        import tarfile
        with tarfile.open(packfile_name, "w:gz") as tar_ref:
            for ifile in packfile_list:
                tar_ref.add(ifile)
    else:
        print("Pack type \'%s\' is not supported!" % pack_type)
        return None
    return packfile_name