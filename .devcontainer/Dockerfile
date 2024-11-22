FROM ubuntu:noble

# Update system
RUN apt-get update -y && apt-get upgrade -y

# Install dependencies for Sparselizard
RUN apt-get install -y \
    clang clang-format cmake zsh git curl \
    libopenblas-dev \
    libopenmpi-dev \
    libmumps-dev \
    libmetis-dev \
    petsc-dev \
    slepc-dev \
    libgmsh-dev \
    libomp-dev

# Fix git issues with the container running as root
RUN git config --global --add safe.directory '*'

# Install Oh My Zsh for convenience
RUN chsh -s $(which zsh)
RUN sh -c "$(curl -fsSL https://raw.githubusercontent.com/ohmyzsh/ohmyzsh/master/tools/install.sh)" || true

# Set the default shell to Zsh
CMD ["zsh"]