###########
# BUILDER #
###########
FROM clinicalgenomics/python3.8-venv:1.0 AS builder

ENV PATH="/venv/bin:$PATH"

WORKDIR /app

# Install dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

#########
# FINAL #
#########
FROM python:3.8-slim

LABEL about.home="https://github.com/J35P312/SVDB"
LABEL about.license="MIT License (MIT)"

## Install base dependencies
RUN apt-get update && \
     apt-get -y upgrade && \
     apt-get -y install -y --no-install-recommends build-essential && \
     apt-get clean && \
     rm -rf /var/lib/apt/lists/*

# Do not upgrade to the latest pip version to ensure more reproducible builds
ENV PIP_DISABLE_PIP_VERSION_CHECK=1
ENV PATH="/venv/bin:$PATH"
RUN echo export PATH="/venv/bin:\$PATH" > /etc/profile.d/venv.sh

## Create a non-root user to run commands
RUN groupadd --gid 1000 worker && useradd -g worker --uid 1000 --shell /usr/sbin/nologin --create-home worker

# Copy virtual environment from builder
COPY --chown=worker:worker --from=builder /venv /venv

WORKDIR /home/worker/app
COPY --chown=worker:worker . /home/worker/app

# Install only Scout app
RUN pip install --no-cache-dir -e .

# Run the app as non-root user
USER worker

ENTRYPOINT [ "svdb" ]
