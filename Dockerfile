###########
# BUILDER #
###########
FROM python:3.12-slim AS builder

RUN python -m venv /venv
ENV PATH="/venv/bin:$PATH"

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

#########
# FINAL #
#########
FROM python:3.12-slim

LABEL about.home="https://github.com/J35P312/SVDB"
LABEL about.license="MIT License (MIT)"

# Upgrade OS packages to reduce known CVEs in the base image, then install
# C build tools needed for optional Cython compilation
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends build-essential && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV PIP_DISABLE_PIP_VERSION_CHECK=1
ENV PATH="/venv/bin:$PATH"
RUN echo export PATH="/venv/bin:\$PATH" > /etc/profile.d/venv.sh

# Non-root user
RUN groupadd --gid 1000 worker && \
    useradd -g worker --uid 1000 --shell /usr/sbin/nologin --create-home worker

# Copy virtual environment from builder stage
COPY --chown=worker:worker --from=builder /venv /venv

WORKDIR /home/worker/app
COPY --chown=worker:worker . /home/worker/app

# Install the package (non-editable for production)
RUN pip install --no-cache-dir .

USER worker

ENTRYPOINT ["svdb"]
