# `rscolorq` changelog

## Version 0.1.2 - 2021-02-20
Previous behavior relied on wrapping overflow in some places which could panic
in debug mode. This behavior has been changed to use checked arithmetic
operations. The algorithm results will differ slightly from the previous
version.

[#6][6] -  Add checked usize arithmetic for Matrix indexing, make clippy fixes  
[#5][5] -  Add safe wrapping usize arithmetic for Matrix2d indexing

## Version 0.1.1 - 2020-11-16
Upgrade the dependencies because an upstream crate made it impossible to compile
the crate and install it from crates.io using `cargo`.

## Version 0.1.0 - 2020-10-20
- Initial Commit

[6]: https://github.com/okaneco/rscolorq/pull/6
[5]: https://github.com/okaneco/rscolorq/pull/5
