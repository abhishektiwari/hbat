"""
File upload component for HBAT web interface.

This module provides a UI component for uploading PDB files.
"""

from typing import Callable, Optional

from nicegui import ui


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

                # Check file size (1 MB = 1048576 bytes)
                MAX_FILE_SIZE = 1 * 1024 * 1024  # 1 MB in bytes
                file_size = len(content)

                if file_size > MAX_FILE_SIZE:
                    size_mb = file_size / (1024 * 1024)
                    ui.notify(
                        f"File too large: {size_mb:.2f} MB. Maximum allowed: 1 MB",
                        type="negative",
                        position="top-left"
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
                ui.notify(f"Upload failed: {str(error)}", type="negative", position="top-left")
                print(f"Upload error: {error}")
                import traceback
                traceback.print_exc()

        async def download_pdb():
            """Download PDB file from RCSB using PDB ID."""
            pdb_id = self.pdb_id_input.value.strip().lower()
            if not pdb_id:
                ui.notify("Please enter a PDB ID", type="warning", position="top-left")
                return

            if len(pdb_id) != 4:
                ui.notify("PDB ID must be 4 characters", type="warning", position="top-left")
                return

            # Validate PDB ID format (alphanumeric only, security check)
            if not pdb_id.isalnum():
                ui.notify("PDB ID must contain only letters and numbers", type="warning", position="top-left")
                return

            try:
                import urllib.request
                import urllib.parse

                # Safely construct URL with validation
                safe_pdb_id = urllib.parse.quote(pdb_id, safe='')
                url = f"https://files.rcsb.org/download/{safe_pdb_id}.pdb"
                ui.notify(f"Downloading {pdb_id} from RCSB PDB...", type="info", position="top-left")

                # Add timeout for security (30 seconds)
                with urllib.request.urlopen(url, timeout=30) as response:  # nosec B310
                    content = response.read()

                # Check file size (1 MB = 1048576 bytes)
                MAX_FILE_SIZE = 1 * 1024 * 1024  # 1 MB in bytes
                file_size = len(content)

                if file_size > MAX_FILE_SIZE:
                    size_mb = file_size / (1024 * 1024)
                    ui.notify(
                        f"Downloaded file too large: {size_mb:.2f} MB. Maximum allowed: 1 MB",
                        type="negative",
                        position="top-left"
                    )
                    return

                filename = f"{pdb_id}.pdb"

                if self.upload_label:
                    self.upload_label.text = f"✓ Downloaded: {filename}"
                    self.upload_label.classes(replace="text-positive")

                self.file_uploaded = True

                if self.on_file_upload:
                    self.on_file_upload(filename, content)

                ui.notify(f"Downloaded {filename} from RCSB PDB", type="positive", position="top-left")
            except Exception as error:
                ui.notify(f"Download failed: {str(error)}", type="negative", position="top-left")
                print(f"Download error: {error}")

        # Store download function for external use
        self._download_pdb_func = download_pdb

        # Option 1: Upload file
        with ui.card().classes("w-full q-pa-md"):
            ui.label("Option 1: Upload PDB File").classes("text-h6")
            ui.upload(
                label="Choose PDB File",
                on_upload=handle_upload,
                auto_upload=True,
            ).props('accept=".pdb"').classes("w-full")

        ui.label("OR").classes("text-center text-bold q-my-md")

        # Option 2: Download from PDB
        with ui.card().classes("w-full q-pa-md"):
            ui.label("Option 2: Download from RCSB PDB").classes("text-h6")
            with ui.row().classes("w-full items-center gap-2"):
                self.pdb_id_input = ui.input(
                    label="PDB ID",
                    placeholder="e.g., 1BHL",
                    validation={"Must be 4 characters": lambda v: len(v.strip()) == 4 if v else True},
                ).props("maxlength=4").classes("flex-1")
                ui.button(
                    "Download",
                    icon="download",
                    on_click=download_pdb,
                ).props("color=primary")

            # Example PDB IDs
            ui.label("Example PDB IDs:").classes("text-caption q-mt-sm")
            with ui.row().classes("gap-1 q-mt-xs flex-wrap"):
                example_ids = ["6RSA", "1CRN", "4HHB", "2PTC", "2IZF", "4LAZ"]
                for pdb_id in example_ids:
                    ui.button(
                        pdb_id,
                        on_click=lambda pid=pdb_id: self.pdb_id_input.set_value(pid)
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
**Supported Files:**

- **PDB files** (.pdb)

- Protein Data Bank format
- Files can be downloaded from [RCSB PDB](https://www.rcsb.org/)
- Maximum file size: 1 MB

**PDB IDs:**

- Enter a 4-character PDB ID (e.g., 1BHL, 1GAI, 1UBI)
- Files are downloaded directly from RCSB PDB
- Common examples: 1BHL (Hemoglobin), 1UBI (Ubiquitin)
            """
            )

    def reset(self):
        """Reset the upload panel to initial state."""
        if self.pdb_id_input:
            self.pdb_id_input.value = ""
        if self.upload_label:
            self.upload_label.text = "No file loaded"
            self.upload_label.classes(replace="text-caption text-grey q-mt-md")
        self.file_uploaded = False
