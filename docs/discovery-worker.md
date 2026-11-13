# Deterministic discovery worker

`discovery_worker.py tick --config CONFIG.json` runs the existing Planet Finder workflow with at most four TICs, three sectors per TIC, 48 MiB of new FITS and seven minutes of subprocess time. Scheduled runs use the Pakistan calendar date. Use `--order ORDER-ID` for a separate bounded work order. An order is a config and a stable ID. Repeating a completed slot or order returns its saved receipt. An interrupted slot reserves its transfer allowance and needs inspection before any new order. The workflow retries one MAST connection failure or timeout after two seconds. FITS and measurement caches can be reused across runs.

On Windows, use `install_discovery_task.ps1 -Python PATH -Config CONFIG.json` under your normal user account. The machine needs Pakistan Standard Time, an interactive user session, a working Python environment and network access at 16:00 Tuesday through Sunday. In this session the sandbox account could not register the task. GitHub Actions also supports manual runs and the same UTC schedule. Its scheduled jobs are gated on `PLANET_FINDER_DISCOVERY=enabled`. Hosted cache and artifact retention are limited, so preserve the queue in durable storage before relying on hosted cache for dispatch deduplication.

## Review rules

A daily actionable candidate needs a blind measured period, at least three supported cycles, zero sensitive missing cycles, complete ingestion and catalog checks, and no prior match or substantive vetting flags. This gate queues a review. It cannot establish novelty, and the default community catalog has a known cutoff. Ambiguous, inconclusive and otherwise flagged signals with at least three supporting cycles stay available for Monday uncertainty review. All measurements remain in the underlying report. The same computation failure in two distinct reports becomes actionable. Archive transport, missing observations and ordinary negative results stay in run receipts.

- `python discovery_worker.py claim` leases up to three actionable packets per Pakistan day.
- `python discovery_worker.py claim --monday-batch` leases up to three uncertainty packets during the existing Monday session.
- `python discovery_worker.py ack EVENT_ID "review outcome"` acknowledges a leased packet. An expired claim becomes uncertain delivery. Reconcile its outcome before another attempt. If you verified no handoff occurred, add `--never-delivered` for one retry after five minutes.

The queue uses SQLite transactions, evidence fingerprints and a one-day cooldown per target after acknowledgment. Packet bytes are a context-size proxy. The worker issues zero model requests. A reviewer can request a different TIC or bounds through an explicit work order. Cross-survey work stays with the reviewer because the worker currently runs the TESS pipeline.

The Codex app's existing heartbeat launches a model each time. The available app thread tool can be called from an agent session, and no worker-callable authenticated dispatch was verified. The installed CLI failed to load its account configuration under the worker identity. Both existing heartbeats therefore remain in place. The GitHub schedule is opt-in so this implementation can be evaluated without creating competing daily work.
