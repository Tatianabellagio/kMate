# `docs/wiki/`: source for the GitHub wiki

These pages are the **user-facing front door**: install, inputs, build a panel, run,
read the output, troubleshoot. They are kept here, in the repo, so they are versioned
and reviewed with the code. The wiki is a published copy.

Deeper technical material stays in `docs/` and is linked from the wiki rather than
duplicated: `ALGORITHM.md` (the math), `PIPELINE_STATE.md` (production recipe),
`BUILDING_A_PANEL.md` (the long-form version of the wiki page).

## Publishing

A GitHub wiki only becomes a git repository **after its first page is created in the
web UI**. One-time setup:

1. Go to <https://github.com/Tatianabellagio/kMate/wiki> and click **Create the first
   page** (any content; it gets overwritten).
2. Then push these pages:

```bash
git clone https://github.com/Tatianabellagio/kMate.wiki.git /tmp/kmate.wiki
cp docs/wiki/*.md /tmp/kmate.wiki/
cd /tmp/kmate.wiki && rm -f README.md          # not a wiki page
git add -A && git commit -m "Publish kMate wiki" && git push
```

`_Sidebar.md` becomes the navigation on every page. File names are page names:
`Building-a-panel.md` → the **Building a panel** page, linked as `[...](Building-a-panel)`.

## Keeping it honest

When a flag, default or requirement changes, update the wiki page in the same commit
as the code. A wiki that drifts from the program is worse than no wiki, because
readers cannot tell it has.
