from pathlib import Path
from datetime import datetime
import json
import re
import shutil
import time
import uuid

import requests
import streamlit as st


# ==================================================
# Public Test Version
# ==================================================

APP_VERSION = "Public Test v1.0"

MAX_UPLOAD_MB = 50
MAX_UPLOAD_BYTES = MAX_UPLOAD_MB * 1024 * 1024

BASE_SESSION_FOLDER = Path("public_test_sessions")
BASE_SESSION_FOLDER.mkdir(exist_ok=True)

WATER_NAMES = {"HOH", "WAT", "H2O", "TIP3", "SOL"}
ION_NAMES = {
    "NA", "CL", "K", "CA", "MG", "ZN", "MN", "FE",
    "CU", "CO", "NI", "CD", "HG", "BR", "IOD", "LI",
    "RB", "SR", "BA"
}


SYNTHETIC_PDB = """HEADER    SYNTHETIC PROTEIN-LIGAND COMPLEX
TITLE     SAMPLE PDB FOR CHARMM-GUI STREAMLINER PUBLIC TESTING
ATOM      1  N   ALA A   1      11.104  13.207   9.601  1.00 20.00           N
ATOM      2  CA  ALA A   1      12.560  13.350   9.700  1.00 20.00           C
ATOM      3  C   ALA A   1      13.140  12.030  10.200  1.00 20.00           C
ATOM      4  O   ALA A   1      12.540  11.000  10.100  1.00 20.00           O
ATOM      5  CB  ALA A   1      13.100  14.500  10.500  1.00 20.00           C
ATOM      6  N   GLY A   2      14.300  12.060  10.760  1.00 20.00           N
ATOM      7  CA  GLY A   2      14.980  10.850  11.260  1.00 20.00           C
ATOM      8  C   GLY A   2      16.450  11.120  11.600  1.00 20.00           C
ATOM      9  O   GLY A   2      17.240  10.200  11.700  1.00 20.00           O
HETATM   10  C1  DOP B   1      18.200  10.900  12.200  1.00 20.00           C
HETATM   11  C2  DOP B   1      19.300  11.200  12.800  1.00 20.00           C
HETATM   12  O1  DOP B   1      20.100  10.300  13.100  1.00 20.00           O
HETATM   13 NA   NA  C   1      15.500  15.100  10.500  1.00 20.00          NA
HETATM   14  O   HOH D   1      10.500  15.600   8.900  1.00 20.00           O
END
"""


# ==================================================
# Session and safety helpers
# ==================================================

def get_session_id():
    if "session_id" not in st.session_state:
        st.session_state["session_id"] = str(uuid.uuid4())[:8]
    return st.session_state["session_id"]


def get_session_folder():
    session_id = get_session_id()
    folder = BASE_SESSION_FOLDER / session_id
    folder.mkdir(parents=True, exist_ok=True)
    return folder


def safe_filename(filename):
    filename = filename.strip().replace(" ", "_")
    filename = re.sub(r"[^A-Za-z0-9_.-]", "", filename)

    if not filename.lower().endswith(".pdb"):
        filename = filename + ".pdb"

    return filename


def clear_current_session_files():
    folder = get_session_folder()
    if folder.exists():
        shutil.rmtree(folder)

    folder.mkdir(parents=True, exist_ok=True)

    keys_to_clear = [
        "pdb_file",
        "cleaned_file",
        "job_id",
        "input_type",
        "original_filename",
        "job_history",
    ]

    for key in keys_to_clear:
        if key in st.session_state:
            del st.session_state[key]


# ==================================================
# PDB functions
# ==================================================

def get_residue_name(line):
    if len(line) >= 20:
        residue = line[17:20].strip().upper()
        if residue:
            return residue

    parts = line.split()
    if len(parts) >= 4:
        return parts[3].upper()

    return ""


def create_synthetic_pdb():
    folder = get_session_folder()
    file_path = folder / "synthetic_protein_ligand.pdb"
    file_path.write_text(SYNTHETIC_PDB, encoding="utf-8")

    st.session_state["original_filename"] = "synthetic_protein_ligand.pdb"
    return file_path


def save_uploaded_pdb(uploaded_file):
    if uploaded_file.size > MAX_UPLOAD_BYTES:
        raise ValueError(f"File is too large. Maximum allowed size is {MAX_UPLOAD_MB} MB.")

    filename = safe_filename(uploaded_file.name)

    folder = get_session_folder()
    file_path = folder / filename

    with open(file_path, "wb") as file:
        file.write(uploaded_file.getbuffer())

    st.session_state["original_filename"] = filename
    return file_path


def fetch_pdb_from_rcsb(pdb_id):
    pdb_id = pdb_id.strip().upper()

    if not re.fullmatch(r"[A-Za-z0-9]{4}", pdb_id):
        raise ValueError("Please enter a valid 4-character PDB ID, for example 2HAC.")

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url, timeout=30)

    if response.status_code != 200:
        raise ValueError(f"Could not fetch PDB ID {pdb_id}. Please check the ID.")

    folder = get_session_folder()
    file_path = folder / f"{pdb_id}.pdb"
    file_path.write_text(response.text, encoding="utf-8")

    st.session_state["original_filename"] = f"{pdb_id}.pdb"
    return file_path


def validate_pdb_basic(file_path):
    atom_or_hetatm_found = False

    with open(file_path, "r", encoding="utf-8", errors="ignore") as file:
        for line in file:
            record = line[:6].strip().upper()
            if record in {"ATOM", "HETATM"}:
                atom_or_hetatm_found = True
                break

    if not atom_or_hetatm_found:
        raise ValueError("This file does not look like a valid PDB file. No ATOM or HETATM records were found.")


def analyse_pdb(file_path):
    validate_pdb_basic(file_path)

    protein_atoms = 0
    hetatm_lines = 0
    ligands = set()
    ions = set()
    waters = set()
    chains = set()
    residues = set()

    with open(file_path, "r", encoding="utf-8", errors="ignore") as file:
        for line in file:
            record = line[:6].strip().upper()

            if record == "ATOM":
                protein_atoms += 1

                if len(line) > 22:
                    chain_id = line[21].strip()
                    if chain_id:
                        chains.add(chain_id)

                residue = get_residue_name(line)
                if residue:
                    residues.add(residue)

            elif record == "HETATM":
                hetatm_lines += 1
                residue_name = get_residue_name(line)

                if residue_name in WATER_NAMES:
                    waters.add(residue_name)
                elif residue_name in ION_NAMES:
                    ions.add(residue_name)
                else:
                    ligands.add(residue_name)

    return {
        "protein_atoms": protein_atoms,
        "hetatm_lines": hetatm_lines,
        "chains": sorted(chains),
        "residue_types": sorted(residues),
        "ligands": sorted(ligands),
        "ions": sorted(ions),
        "waters": sorted(waters),
        "has_protein": protein_atoms > 0,
        "has_ligand": len(ligands) > 0,
        "has_ions": len(ions) > 0,
        "has_water": len(waters) > 0,
    }


def clean_pdb(input_path, keep_protein, keep_ligand, keep_ions, keep_water):
    folder = get_session_folder()

    original_stem = Path(st.session_state.get("original_filename", "input.pdb")).stem
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    output_path = folder / f"{original_stem}_cleaned_{timestamp}.pdb"

    cleaned_lines = []

    with open(input_path, "r", encoding="utf-8", errors="ignore") as file:
        for line in file:
            record = line[:6].strip().upper()

            if record == "ATOM":
                if keep_protein:
                    cleaned_lines.append(line)

            elif record == "HETATM":
                residue_name = get_residue_name(line)

                if residue_name in WATER_NAMES:
                    if keep_water:
                        cleaned_lines.append(line)

                elif residue_name in ION_NAMES:
                    if keep_ions:
                        cleaned_lines.append(line)

                else:
                    if keep_ligand:
                        cleaned_lines.append(line)

    cleaned_lines.append("END\n")

    if len(cleaned_lines) <= 1:
        raise ValueError("The cleaned PDB would be empty. Please keep at least one valid component.")

    output_path.write_text("".join(cleaned_lines), encoding="utf-8")
    return output_path


def preview_file(file_path, max_lines=25):
    with open(file_path, "r", encoding="utf-8", errors="ignore") as file:
        lines = file.readlines()[:max_lines]
    return "".join(lines)


# ==================================================
# Demo job functions
# ==================================================

def get_job_history():
    if "job_history" not in st.session_state:
        st.session_state["job_history"] = []
    return st.session_state["job_history"]


def demo_submit_job():
    return f"DEMO-{int(time.time())}"


def add_job_to_session_history(job_record):
    history = get_job_history()
    history.append(job_record)
    st.session_state["job_history"] = history


# ==================================================
# Streamlit UI
# ==================================================

st.set_page_config(
    page_title="CHARMM-GUI Streamliner Public Test",
    layout="wide"
)

st.title("CHARMM-GUI Streamliner — Public Test")

st.caption(APP_VERSION)

st.info(
    "This public test version accepts real PDB uploads, analyses protein/ligand/ions/water, "
    "creates a cleaned PDB, and lets users download the cleaned file. "
    "It does not submit jobs to CHARMM-GUI yet."
)

st.warning(
    "Please do not upload confidential or sensitive structures during public testing. "
    "Uploaded files are used only for this session, but this prototype is not a secure production system."
)


# ==================================================
# Sidebar
# ==================================================

st.sidebar.header("Testing Options")

st.sidebar.success("Demo mode is always ON for public testing.")

session_id = get_session_id()
st.sidebar.write(f"Session ID: `{session_id}`")

if st.sidebar.button("Clear my session files", key="clear_session_button"):
    clear_current_session_files()
    st.sidebar.success("Session files cleared.")
    st.rerun()


# ==================================================
# Step 1: Input
# ==================================================

st.header("Step 1: Choose PDB input")

input_choice = st.radio(
    "Choose input type:",
    ["Use synthetic test PDB", "Upload real PDB file", "Fetch PDB from RCSB ID"],
    key="input_choice_radio"
)


if input_choice == "Use synthetic test PDB":
    if st.button("Create synthetic PDB file", key="create_synthetic_button"):
        try:
            pdb_file = create_synthetic_pdb()
            st.session_state["pdb_file"] = str(pdb_file)
            st.session_state["input_type"] = "Synthetic test PDB"

            st.success("Synthetic PDB file created.")
            st.write(f"File location: `{pdb_file.name}`")
        except Exception as error:
            st.error(str(error))


elif input_choice == "Upload real PDB file":
    uploaded_file = st.file_uploader(
        "Upload a .pdb file",
        type=["pdb"],
        key="real_pdb_uploader",
        help=f"Maximum recommended test size: {MAX_UPLOAD_MB} MB."
    )

    if uploaded_file is not None:
        st.write(f"Selected file: `{uploaded_file.name}`")
        st.write(f"File size: `{uploaded_file.size / 1024 / 1024:.2f} MB`")

        if st.button("Save uploaded PDB file", key="save_uploaded_button"):
            try:
                pdb_file = save_uploaded_pdb(uploaded_file)
                validate_pdb_basic(pdb_file)

                st.session_state["pdb_file"] = str(pdb_file)
                st.session_state["input_type"] = "Uploaded real PDB"

                st.success("Uploaded PDB file saved.")
                st.write(f"File location: `{pdb_file.name}`")
            except Exception as error:
                st.error(str(error))


elif input_choice == "Fetch PDB from RCSB ID":
    pdb_id = st.text_input(
        "Enter 4-character PDB ID",
        placeholder="Example: 2HAC",
        max_chars=4,
        key="rcsb_pdb_id"
    )

    if st.button("Fetch PDB", key="fetch_pdb_button"):
        try:
            pdb_file = fetch_pdb_from_rcsb(pdb_id)
            validate_pdb_basic(pdb_file)

            st.session_state["pdb_file"] = str(pdb_file)
            st.session_state["input_type"] = "Fetched from RCSB"

            st.success(f"PDB {pdb_id.upper()} fetched successfully.")
            st.write(f"File location: `{pdb_file.name}`")
        except Exception as error:
            st.error(str(error))


# ==================================================
# Step 2: Analyse
# ==================================================

st.header("Step 2: Analyse PDB file")

if "pdb_file" not in st.session_state:
    st.warning("Please create, upload, or fetch a PDB file first.")

else:
    try:
        pdb_path = st.session_state["pdb_file"]
        summary = analyse_pdb(pdb_path)

        st.info(f"Current PDB file: `{Path(pdb_path).name}`")

        col1, col2, col3, col4, col5 = st.columns(5)

        col1.metric("Protein atoms", summary["protein_atoms"])
        col2.metric("HETATM lines", summary["hetatm_lines"])
        col3.metric("Chains", len(summary["chains"]))
        col4.metric("Ligands", len(summary["ligands"]))
        col5.metric("Water types", len(summary["waters"]))

        st.subheader("Detected content")

        st.write("Chains found:", summary["chains"])
        st.write("Ligands found:", summary["ligands"])
        st.write("Ions found:", summary["ions"])
        st.write("Waters found:", summary["waters"])

        with st.expander("Preview first 25 lines of PDB file"):
            st.code(preview_file(pdb_path), language="text")

    except Exception as error:
        st.error(str(error))


# ==================================================
# Step 3: Clean
# ==================================================

st.header("Step 3: Choose what to keep")

if "pdb_file" not in st.session_state:
    st.warning("Create, upload, or fetch a PDB file first.")

else:
    try:
        summary = analyse_pdb(st.session_state["pdb_file"])

        keep_protein = st.checkbox(
            "Keep protein / amino acids",
            value=summary["has_protein"],
            key="keep_protein"
        )

        keep_ligand = st.checkbox(
            "Keep ligand / drug",
            value=summary["has_ligand"],
            key="keep_ligand"
        )

        keep_ions = st.checkbox(
            "Keep ions",
            value=summary["has_ions"],
            key="keep_ions"
        )

        keep_water = st.checkbox(
            "Keep water",
            value=summary["has_water"],
            key="keep_water"
        )

        if st.button("Create cleaned PDB", key="create_cleaned_button"):
            try:
                cleaned_file = clean_pdb(
                    input_path=st.session_state["pdb_file"],
                    keep_protein=keep_protein,
                    keep_ligand=keep_ligand,
                    keep_ions=keep_ions,
                    keep_water=keep_water,
                )

                st.session_state["cleaned_file"] = str(cleaned_file)

                st.success("Cleaned PDB file created.")
                st.write(f"Cleaned file: `{cleaned_file.name}`")

                with open(cleaned_file, "rb") as file:
                    st.download_button(
                        label="Download cleaned PDB",
                        data=file,
                        file_name=cleaned_file.name,
                        mime="chemical/x-pdb",
                        key="download_cleaned_button"
                    )

            except Exception as error:
                st.error(str(error))

    except Exception as error:
        st.error(str(error))


# ==================================================
# Step 4: Demo job
# ==================================================

st.header("Step 4: Demo job submission")

if "cleaned_file" not in st.session_state:
    st.warning("Create the cleaned PDB first.")

else:
    if st.button("Submit demo job", key="submit_demo_button"):
        job_id = demo_submit_job()
        st.session_state["job_id"] = job_id

        job_record = {
            "job_id": job_id,
            "status": "DONE",
            "mode": "Demo",
            "input_type": st.session_state.get("input_type", "Unknown"),
            "input_file": Path(st.session_state.get("pdb_file", "Unknown")).name,
            "cleaned_file": Path(st.session_state.get("cleaned_file", "Unknown")).name,
            "created_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }

        add_job_to_session_history(job_record)

        st.success("Demo job submitted.")
        st.write(f"Job ID: `{job_id}`")


# ==================================================
# Step 5: Status
# ==================================================

st.header("Step 5: Check job status")

if "job_id" not in st.session_state:
    st.warning("Submit a demo job first.")

else:
    if st.button("Check job status", key="check_status_button"):
        st.success("Job status: DONE")
        st.write("This is demo status. Real CHARMM-GUI job submission will be added after API confirmation.")


# ==================================================
# Step 6: Session history
# ==================================================

st.header("Step 6: Session job history")

history = get_job_history()

if not history:
    st.info("No jobs in this session yet.")

else:
    st.dataframe(history, use_container_width=True)

    history_json = json.dumps(history, indent=4)

    st.download_button(
        label="Download session job history JSON",
        data=history_json,
        file_name="session_job_history.json",
        mime="application/json",
        key="download_history_button"
    )


# ==================================================
# Footer
# ==================================================

st.divider()

st.subheader("Current public-test scope")

st.write(
    """
This version is ready for external testers to try with real `.pdb` files.

It supports:
- synthetic PDB testing
- real PDB upload
- RCSB PDB fetch by ID
- detection of protein atoms, chains, ligands, ions, and water
- cleaned PDB generation
- cleaned PDB download
- session-only demo job history

It does not yet:
- submit jobs to real CHARMM-GUI
- store CHARMM-GUI credentials
- perform ligand parameterisation
- run molecular dynamics
"""
)