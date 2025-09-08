# Use a standard Python image as a base
FROM python:3.8-slim

# Set the working directory in the container
WORKDIR /app

# Copy the local repository content to the container
COPY . .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Install ExonSurfer
RUN python setup.py install

# Set the entrypoint to the exon_surfer.py script
ENTRYPOINT ["exon_surfer.py"]
