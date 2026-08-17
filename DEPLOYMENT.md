# Deploying CreMAP

## 1. Publish the source on GitHub

The repository contains application source and an `renv.lock` dependency lock.
Local R libraries, credentials, histories, caches, and `*.rds` datasets are
excluded from Git.

Create an empty public GitHub repository, then run from this directory:

```bash
git add .
git commit -m "Prepare CreMAP for public deployment"
git remote add origin https://github.com/YOUR_USER/YOUR_REPO.git
git push -u origin main
```

Do not commit the shinyapps.io token or secret. `rsconnect::setAccountInfo()`
stores those in the user account configuration outside this repository.

## 2. Register a shinyapps.io account locally

Install the deployment client outside the app's runtime library:

```bash
R --vanilla -e 'install.packages("rsconnect")'
```

In the shinyapps.io dashboard, open **Account > Tokens > Show > Show Secret**,
then run the `setAccountInfo(...)` command shown there. Do not save that command
in a repository file or shell script.

## 3. Optional: include real reference atlases

Without the two reference files, the deployment intentionally uses synthetic
mouse and human demonstration data. To deploy the real references, build them
locally first:

```bash
Rscript scripts/install_deps.R
Rscript scripts/download_reference_data.R
```

The resulting `reference/*.rds` files remain ignored by Git, but the deployment
bundle includes them. Check their combined size before deploying. shinyapps.io
bundle limits depend on the subscription plan.

## 4. Deploy

Optionally select the account and application name, then run the deployment
script with vanilla R:

```bash
export SHINYAPPS_ACCOUNT="YOUR_ACCOUNT"
export SHINYAPPS_APP_NAME="cremap"
R --vanilla -f scripts/deploy_shinyapps.R
```

The application will be available at:

```text
https://YOUR_ACCOUNT.shinyapps.io/cremap/
```

The `.rscignore` file keeps development-only files out of the application
bundle. `renv.lock` supplies the hosted R package versions.

## Hosted defaults

- Maximum compressed upload: 100 MB (`CREMAP_MAX_UPLOAD_MB`).
- Maximum expanded Cell Ranger ZIP: 500 MB (`CREMAP_MAX_UNCOMPRESSED_MB`).
- Live ExperimentHub atlas downloads are disabled automatically on
  shinyapps.io. Bundle local RDS files for real defaults.
- MouseMine requests remain live and require outbound HTTPS.
- Uploaded data and manual session changes are ephemeral.

The upload and analysis features can consume substantial memory. Start by
testing synthetic data on the default instance. A 2–4 GB instance is more
appropriate for public Seurat/Cell Ranger uploads.

## 5. Embed in Google Sites

In Google Sites choose **Insert > Embed > By URL** and use the shinyapps.io URL.
A full-page embed works best for CreMAP's navigation and plots. Also provide a
normal link that opens the application in a new tab.
