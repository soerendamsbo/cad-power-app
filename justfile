# justfile — task runner for WP2 power analysis app
# Usage: just <recipe>
# List all recipes: just --list

app_dir := "shiny"

# Launch the Shiny app (requires results.rds — run `just precompute` first)
app:
    Rscript -e "shiny::runApp('{{app_dir}}/app.R', launch.browser = TRUE)"

# Run the precompute grid (writes shiny/results.rds, ~15–25 min)
precompute:
    Rscript {{app_dir}}/precompute.R

# Precompute then launch
all: precompute app
