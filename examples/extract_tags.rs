// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::{Tagger, TaggerOptions};
use openms::comparison::Tolerance;

fn main() -> openms::Result<()> {
    // Rounded source P/E/P/T increments, independent of the computed mass table.
    let mut positions = vec![150.0];
    for increment in [97.0527, 129.0426, 97.0527, 101.0477] {
        positions.push(positions.last().unwrap() + increment);
    }
    let mut options = TaggerOptions::new(2, Tolerance::Absolute(0.02));
    options.max_tag_length = 4;
    let tags = Tagger::new(options)?.get_tags(&positions)?;
    for tag in &tags {
        println!("{tag}");
    }
    println!(
        "{} distinct tags from {} measured positions",
        tags.len(),
        positions.len()
    );
    Ok(())
}
