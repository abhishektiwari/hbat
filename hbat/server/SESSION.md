# HBAT Web Server - Session Management

## Overview

The HBAT web server now uses UUID-based session management to isolate data for multiple concurrent users. Each user gets their own dedicated session directory where all uploaded files and analysis results are stored.

## Session Lifecycle

### Creation
- A new session is automatically created when a user uploads their first PDB file
- Each session gets a unique UUID (e.g., `a1b2c3d4-e5f6-7890-abcd-ef1234567890`)
- A `.session_timestamp` file tracks when the session was created

### Active Period
- Sessions remain valid for 24 hours from creation
- All user files are stored in their session directory: `uploads/sessions/<uuid>/`

### Cleanup
- On startup: All expired sessions (>24 hours old) are removed
- Periodic: Every 6 hours, expired sessions are automatically cleaned up
- Manual: When a user clicks "Start Over", their session is immediately deleted

## Configuration

### Session Timeout
Default timeout is 24 hours. To change it, modify `hbat/server/app.py`:

```python
session_manager = SessionManager(SESSIONS_BASE_DIR, session_timeout_hours=48)  # 48 hours
```

### Cleanup Interval
Default cleanup runs every 6 hours. To change it, modify the `periodic_cleanup()` function in `hbat/server/app.py`:

```python
await asyncio.sleep(12 * 3600)  # 12 hours instead of 6
```

## Docker Compose Example

```yaml
version: '3.8'

services:
  hbat-web:
    image: hbat:latest
    volumes:
      - ./uploads:/app/uploads  # Single volume mount for all user data
    ports:
      - "8080:8080"
    environment:
      - HBAT_ENV=production
      - HBAT_HOST=0.0.0.0
      - HBAT_PORT=8080
```

## Session Files

Each session may contain:
- Uploaded PDB files: Original files uploaded by the user
- Fixed PDB files: PDB files after hydrogen addition (if PDB fixing was applied)
- Result files: JSON, TXT, CSV exports of analysis results
- Graph files: PNG/SVG files for cooperativity chain visualizations
- Timestamp: `.session_timestamp` file marking session creation time

## Security Considerations

1. Random UUIDs: Session IDs are generated using `uuid.uuid4()` for cryptographic randomness
2. Path isolation: Users cannot specify session IDs or access other sessions
3. File size limits: Uploads are limited to 1 MB per file
4. Automatic cleanup: Old data is automatically removed after timeout

## Implementation Details

### SessionManager Class
Located in `hbat/server/session.py`, provides:
- `create_session()`: Generate new UUID-based session
- `get_session_dir(session_id)`: Get path to session directory
- `session_exists(session_id)`: Check if session is valid
- `cleanup_expired_sessions()`: Remove all expired sessions
- `delete_session(session_id)`: Delete specific session
- `get_session_file_path(session_id, filename)`: Get full path for a file

### Integration
The `HBATWebApp` class in `hbat/server/app.py`:
- Creates a session on first file upload
- Stores all files in the session directory
- Passes session directory to results panel for graph generation
- Cleans up session on "Start Over"

## Monitoring

To monitor active sessions:

```bash
# Count active sessions
ls -1 uploads/sessions/ | wc -l

# Show session ages
find uploads/sessions/ -name ".session_timestamp" -exec stat -f "%Sm %N" {} \;

# Manually trigger cleanup (requires Python)
python -c "from hbat.server.session import SessionManager; from pathlib import Path; sm = SessionManager(Path('uploads/sessions')); print(f'Cleaned {sm.cleanup_expired_sessions()} sessions')"
```
