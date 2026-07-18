.PHONY: setup pipeline dashboard clean

# Install all project dependencies.
setup:
	pip install -r requirements.txt

# Run the full pipeline: init DB + load data (Part 1), then all analysis (Parts 2-4).
pipeline:
	python load_data.py
	python analysis.py

# Start the interactive dashboard (local server).
dashboard:
	streamlit run dashboard.py

# Remove generated artifacts.
clean:
	rm -f cell_counts.db
	rm -rf outputs
