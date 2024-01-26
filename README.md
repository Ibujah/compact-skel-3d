# 3D skeletonization

Compact 3d skeletonization in rust, mesh conversion to delaunay. Link to the article [here](https://hal.science/hal-04262568) (prepublication on HAL).

## Instructions

### Installing Rust

This program is writen in Rust. You can install the main tool (cargo) [here](https://www.rust-lang.org/tools/install).

### Build and run

Build instructions:
```
cargo build --release
```

Mesh conversion to Delaunay:
```
cargo run --release --bin soft_todelaunay -- --meshinfile ./ressources/hand.obj --objoutfile ./ressources/hand_del.obj
```

Delaunay mesh skeletonization:
```
cargo run --release --bin soft_sheetskeletonization -- --meshinfile ./ressources/hand_del.obj --epsilon 0.01 --pathout ./hand/
```
