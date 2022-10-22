from typing import Optional

import ase.build
import click
from ase.data.pubchem import pubchem_atoms_conformer_search

from . import utils
from .console import console

available_types = click.Choice(["sdf", "xyz", "pdb"])


@click.group()
def app():
    """Get a molecule from various sources and write it to a file.

    Currently available sources:

        - ASE (by ase.build.molecule)
        - PubChem (by CID, name, or SMILES)
        - RDKit (by arbitrary SMILES)

    See help for each subcommand for more information.
    """
    pass


@app.command()
@click.argument("formula")
@click.option("-o", "--output", default=None, help="Output filename")
def from_formula(formula: str, output: Optional[str] = None):
    """Generate molecule from chemical formula (ASE)"""
    if output is None:
        output = f"{formula}.xyz"
    atoms = ase.build.molecule(formula)
    atoms.write(output)


@app.command()
@click.option("--cid", type=int, default=None, help="PubChem CID")
@click.option("--name", type=str, default=None, help="PubChem name")
@click.option("--smiles", type=str, default=None, help="SMILES string")
@click.option("--show-info-only", is_flag=True, help="Show molecule info and exit.")
@click.option("--all", "all_conformers", is_flag=True, default=False, help="Find all conformers")
@click.option("-o", "--output", default=None, help="Output filename")
@click.option("-f", "--format", "output_format", type=available_types, default="xyz", help="Output format")
def from_pubchem(
    cid: Optional[int] = None,
    name: Optional[str] = None,
    smiles: Optional[str] = None,
    show_info_only: bool = False,
    all_conformers: bool = False,
    output: Optional[str] = None,
    output_format: Optional[str] = "xyz",
):
    """Search molecule by Pubchem API.
    Exactly one of cid, name, or SMILES must be specified.
    """
    namespace, value = utils.extract_field(cid, name, smiles)
    filename, output_format = utils.determine_output(namespace, value, output, output_format)
    compound = utils.pubchem_get_compound(namespace, value)
    utils.show_compound_info(compound)
    if show_info_only:
        return

    if not all_conformers:
        output = f"{namespace}_{value}.xyz" if output is None else output
        # generate from canonical smiles
        mol = utils.mol_from_smiles(compound.canonical_smiles)
        utils.write_mol(mol, filename, output_format)

    else:
        if output_format != "xyz":
            raise ValueError("Currently only xyz output is supported")
        console.print("Finding all optimized conformers...", style="info")
        conformers = pubchem_atoms_conformer_search(cid=compound.cid)
        console.print(f"Found {len(conformers)} conformers", style="info")
        for i, atoms in enumerate(conformers):
            atoms.set_initial_charges(None)
            filename_i = utils.insert_index(filename, i)
            atoms.write(filename_i)


@app.command()
@click.argument("smiles")
@click.option("--optimize", default=True, help="Optimize geometry with force field")
@click.option("-o", "--output", type=str, default=None, help="Output filename")
@click.option(
    "-f",
    "--format",
    "output_format",
    type=click.Choice(["sdf", "xyz", "pdb"]),
    default="xyz",
    help="Output format",
)
def from_smiles(
    smiles: str,
    optimize: bool,
    output: click.File,
    output_format: str,
):
    """Generate molecule from SMILES string (RDKit).
    You can optimize the molecular geometry with a force field,
    but this may be incorrect.
    """
    console.print("Generating molecule from smiles", style="info")
    mol = utils.mol_from_smiles(smiles, optimize)
    filename, output_format = utils.determine_output("smiles", smiles, output, output_format)
    console.print(f"Writing molecule to {filename}", style="info")
    utils.write_mol(mol, filename, output_format)
