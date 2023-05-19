# 3D skeletonization

Compact 3d skeletonization in rust, mesh conversion to delaunay

## Instructions

### Installing Rust

This program is writen in Rust. You can install the main tool (cargo) [here](https://www.rust-lang.org/tools/install).

### Installing dependencies

Linux:
```
sudo apt-get install build essential liblapacke-dev libopenblas-dev
pip install matplotlib
```

MacOS
```
brew install openblas lapack
pip3 install matplotlib
```

### Build and run

Build instructions:
```
cargo build --release
```

Mesh conversion to Delaunay:
```
cargo run --release --bin soft_todelaunay -- --objinfile ./ressources/hand.obj --objoutfile ./ressources/hand_del.obj
```

Delaunay mesh skeletonization:
```
cargo run --release --bin soft_sheetskeletonization -- --objinfile ./ressources/hand_del.obj --epsilon 0.01 --pathout ./hand/
```
