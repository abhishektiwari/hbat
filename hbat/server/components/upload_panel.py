"""
File upload component for HBAT web interface.

This module provides a UI component for uploading PDB files.
"""

from typing import Callable, Optional

from nicegui import ui

from ..utils.validators import validate_pdb_id


class UploadPanel:
    """File upload panel component."""

    def __init__(self, on_file_upload: Optional[Callable] = None):
        """Initialize the upload panel.

        :param on_file_upload: Callback function when file is uploaded
        :type on_file_upload: Optional[Callable]
        """
        self.on_file_upload = on_file_upload
        self.upload_label: Optional[ui.label] = None
        self.pdb_id_input: Optional[ui.input] = None
        self.format_select: Optional[ui.select] = None
        self.upload_widget: Optional[ui.upload] = None
        self.file_uploaded: bool = False
        self._download_pdb_func: Optional[Callable] = None

    async def download_if_pdb_id_provided(self) -> bool:
        """Download PDB file if PDB ID is provided but file not uploaded yet.

        :returns: True if download was attempted and successful, False otherwise
        :rtype: bool
        """
        if self.file_uploaded:
            return True

        if self.pdb_id_input and self.pdb_id_input.value.strip():
            if self._download_pdb_func:
                await self._download_pdb_func()
                return self.file_uploaded

        return False

    def create_ui(self):
        """Create the upload UI components."""

        async def handle_upload(e):
            """Handle file upload event."""
            try:
                # Access file info through e.file
                filename = e.file.name
                # read() is async, so we need to await it
                content = await e.file.read()

                # Check file size (2 MB = 2097152 bytes)
                MAX_FILE_SIZE = 2 * 1024 * 1024  # 2 MB in bytes
                file_size = len(content)

                if file_size > MAX_FILE_SIZE:
                    size_mb = file_size / (1024 * 1024)
                    ui.notify(
                        f"File too large: {size_mb:.2f} MB. Maximum allowed: 2 MB",
                        type="negative",
                        position="top-left",
                    )
                    return

                if self.upload_label:
                    self.upload_label.text = f"✓ Uploaded: {filename}"
                    self.upload_label.classes(replace="text-positive")

                self.file_uploaded = True

                if self.on_file_upload:
                    self.on_file_upload(filename, content)

                ui.notify(f"Uploaded {filename}", type="positive", position="top-left")
            except Exception as error:
                ui.notify(
                    f"Upload failed: {str(error)}", type="negative", position="top-left"
                )
                print(f"Upload error: {error}")
                import traceback

                traceback.print_exc()

        async def download_structure():
            """Download structure file from RCSB using PDB ID."""
            pdb_id = self.pdb_id_input.value.strip().lower()

            # Validate PDB ID format
            is_valid, error_message = validate_pdb_id(pdb_id)
            if not is_valid:
                ui.notify(error_message, type="warning", position="top-left")
                return

            try:
                import urllib.parse
                import urllib.request

                # Get selected format
                file_format = self.format_select.value if self.format_select else "pdb"

                # Safely construct URL with validation
                safe_pdb_id = urllib.parse.quote(pdb_id, safe="")

                if file_format == "pdb":
                    url = f"https://files.rcsb.org/download/{safe_pdb_id}.pdb"
                    filename = f"{pdb_id}.pdb"
                else:  # mmcif
                    url = f"https://files.rcsb.org/download/{safe_pdb_id}.cif"
                    filename = f"{pdb_id}.cif"

                # Validate URL scheme and domain before opening
                parsed = urllib.parse.urlparse(url)
                if parsed.scheme != "https" or not parsed.netloc.endswith("rcsb.org"):
                    ui.notify(
                        "Invalid URL scheme or domain",
                        type="negative",
                        position="top-left",
                    )
                    return

                # Add timeout for security (30 seconds)
                # URL is validated to be https://files.rcsb.org only
                with urllib.request.urlopen(url, timeout=30) as response:  # nosec B310
                    content = response.read()

                # Check file size (2 MB = 2097152 bytes)
                MAX_FILE_SIZE = 2 * 1024 * 1024  # 2 MB in bytes
                file_size = len(content)

                if file_size > MAX_FILE_SIZE:
                    size_mb = file_size / (1024 * 1024)
                    ui.notify(
                        f"Downloaded file too large: {size_mb:.2f} MB. Maximum allowed: 2 MB",
                        type="negative",
                        position="top-left",
                    )
                    return

                if self.upload_label:
                    self.upload_label.text = f"✓ Downloaded: {filename}"
                    self.upload_label.classes(replace="text-positive")

                self.file_uploaded = True

                if self.on_file_upload:
                    self.on_file_upload(filename, content)

                ui.notify(
                    f"Downloaded {filename} from RCSB PDB",
                    type="positive",
                    position="top-left",
                )
            except Exception as error:
                ui.notify(
                    f"Download failed: {str(error)}",
                    type="negative",
                    position="top-left",
                )
                print(f"Download error: {error}")

        # Store download function for external use
        self._download_pdb_func = download_structure

        # Option 1: Upload file
        with ui.card().classes("w-full q-pa-md mt-5"):
            ui.label("Option 1: Upload Structure File").classes("text-h6")
            self.upload_widget = (
                ui.upload(
                    label="Choose PDB or CIF File",
                    on_upload=handle_upload,
                    auto_upload=True,
                )
                .props('accept=".pdb,.cif"')
                .classes("w-full")
            )

        ui.label("OR").classes("text-center text-bold q-my-md")

        # Option 2: Download from PDB
        with ui.card().classes("w-full q-pa-md"):
            ui.label("Option 2: Download from RCSB PDB").classes("text-h6")

            # Format selector
            self.format_select = (
                ui.select(
                    label="Format",
                    value="pdb",
                    options={"pdb": "PDB (.pdb)", "mmcif": "mmCIF (.cif)"},
                )
                .props("outlined")
                .classes("w-full")
            )

            # PDB ID input
            self.pdb_id_input = (
                ui.input(
                    label="PDB ID",
                    placeholder="e.g., 1BHL or pdb_00001abc",
                    validation={
                        "Invalid PDB ID format": lambda v: (
                            validate_pdb_id(v)[0] if v else True
                        )
                    },
                )
                .props("maxlength=12")
                .classes("w-full")
            )

            # Example PDB IDs
            ui.label("Example PDB IDs:").classes("text-caption q-mt-sm")
            with ui.row().classes("gap-1 q-mt-xs flex-wrap"):
                example_ids = ["6RSA", "1CRN", "4HHB", "2PTC", "2IZF", "4LAZ"]
                for pdb_id in example_ids:
                    ui.button(
                        pdb_id,
                        on_click=lambda pid=pdb_id: self.pdb_id_input.set_value(pid),
                    ).props("size=sm outline color=primary")

        # Status label
        self.upload_label = ui.label("No file loaded").classes(
            "text-caption text-grey q-mt-md"
        )

        # Info
        with ui.expansion("Supported Formats & Info", icon="info").classes(
            "w-full q-mt-md"
        ):
            ui.markdown(
                """
**Supported File Formats:**

- **PDB** (.pdb) - Protein Data Bank format
- **mmCIF** (.cif) - Macromolecular Crystallographic Information File

Both formats can be uploaded directly or downloaded from [RCSB PDB](https://www.rcsb.org/)

Maximum file size: 2 MB

**Download from RCSB PDB:**

- Enter a PDB ID using one of these formats:
  - **Traditional**: 4-character ID (e.g., 1BHL, 1GAI, 1UBI) - format: 1 digit + 3 alphanumeric
  - **Extended (wwPDB)**: 12-character ID (e.g., pdb_00001abc) - format: "pdb_" + 8 alphanumeric
- Select desired format (PDB or mmCIF)
- Files are downloaded directly from RCSB PDB
- Common examples: 1BHL (Hemoglobin), 1UBI (Ubiquitin)

**Format Notes:**

- **PDB**: Traditional Protein Data Bank format, widely supported
- **mmCIF**: Modern standard for macromolecular structures, handles larger and more complex structures better
            """
            )

    def reset(self):
        """Reset the upload panel to initial state."""
        # Clear PDB ID input
        if self.pdb_id_input:
            self.pdb_id_input.value = ""

        # Clear upload widget queue/selection using NiceGUI's reset() method
        if self.upload_widget:
            self.upload_widget.reset()

        # Reset upload label
        if self.upload_label:
            self.upload_label.text = "No file loaded"
            self.upload_label.classes(replace="text-caption text-grey q-mt-md")

        # Reset file uploaded flag
        self.file_uploaded = False
