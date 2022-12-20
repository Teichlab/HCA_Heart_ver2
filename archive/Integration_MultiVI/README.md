# Integarting RNA, ATAC, and Multiome data with `MultiVI`

1. SSH to the farm: eg `ssh kk18@farm5-head1`
2. [optional]  Create a new screen session in so you can keep things running even if you loose connection: screen -S multivi. Read about screen here.
3. Run the script `./run_multivi_interactive.sh`
4. You'll see something like
```
[I 11:23:33.327 NotebookApp] Serving notebooks from local directory: /nfs/your/path
[I 11:23:33.327 NotebookApp] Jupyter Notebook 6.4.3 is running at:
[I 11:23:33.327 NotebookApp] http://farm5-gpu0103:1555/?token=...
[I 11:23:33.327 NotebookApp]  or http://127.0.0.1:1555/?token=...
[I 11:23:33.327 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
```
Take note of `http://farm5-gpu0103:1555/` and paste that URL on your preferred web browser. When asked for password use the same one you defined in the script, by default I think is `123`.
5. Make sure you're using the kernel `multiVI container` when running your notebooks
