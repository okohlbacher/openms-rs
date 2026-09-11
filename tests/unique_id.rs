use openms::concept::{HasUniqueId, UniqueId, UniqueIdGenerator};
use std::collections::HashSet;

const SEED: u64 = 546_666_321;
const WORDS: [u64; 6] = [
    4_039_984_684_862_977_299,
    11_561_668_883_169_444_769,
    8_153_960_635_892_418_594,
    12_940_485_248_168_291_983,
    11_522_917_731_873_626_020,
    4_387_255_872_055_054_320,
];

#[test]
fn source_seed_words_and_native_state_copy() {
    let mut random = UniqueIdGenerator::from_seed(SEED);
    assert_eq!(random.seed(), SEED);
    for expected in WORDS {
        assert_eq!(random.get_unique_id(), expected);
    }
    let mut copied = random.clone();
    for _ in 0..1250 {
        assert_eq!(random.get_unique_id(), copied.get_unique_id());
    }
    random.set_seed(SEED);
    for expected in WORDS {
        assert_eq!(random.get_unique_id(), expected);
    }
    random.set_seed(u64::MAX);
    assert_eq!(random.seed(), u64::MAX);
    let mut another = UniqueIdGenerator::from_seed(u64::MAX);
    assert_eq!(random.get_unique_id(), another.get_unique_id());
}

#[test]
fn source_uuid_layout_consumes_two_words_and_keeps_the_shared_sequence() {
    let mut random = UniqueIdGenerator::from_seed(SEED);
    let expected = if cfg!(target_endian = "little") {
        "13b54687-aae8-4038-a183-00c8974c73a0"
    } else {
        "3810e8aa-8746-4513-a073-4c97c80083a1"
    };
    assert_eq!(random.get_uuid(), expected);
    assert_eq!(random.get_unique_id(), WORDS[2]);
    let mut seen = HashSet::new();
    for _ in 0..1000 {
        let uuid = random.get_uuid();
        assert_eq!(uuid.len(), 36);
        for i in [8, 13, 18, 23] {
            assert_eq!(uuid.as_bytes()[i], b'-');
        }
        assert_eq!(uuid.as_bytes()[14], b'4');
        assert!(matches!(uuid.as_bytes()[19], b'8' | b'9' | b'a' | b'b'));
        assert!(uuid.bytes().all(|b| b == b'-' || b.is_ascii_hexdigit()));
        assert!(seen.insert(uuid));
    }
}

#[test]
fn source_interface_values_counts_and_generator_consumption() {
    let mut value = UniqueId::default();
    assert_eq!(value.unique_id(), 0);
    assert_eq!(value.clear_unique_id(), 0);
    assert!(value.has_invalid_unique_id());
    assert!(!UniqueId::is_valid(UniqueId::INVALID));
    assert!(UniqueId::is_valid(1_234_567_890));
    value.set_unique_id(17);
    assert_eq!(value, UniqueId(17));
    assert_eq!(value.clear_unique_id(), 1);
    assert_eq!(value.clear_unique_id(), 0);
    let mut random = UniqueIdGenerator::from_seed(SEED);
    assert_eq!(value.ensure_unique_id(&mut random), 1);
    assert_eq!(value.unique_id(), WORDS[0]);
    assert_eq!(value.ensure_unique_id(&mut random), 0);
    assert_eq!(value.assign_new_unique_id(&mut random), 1);
    assert_eq!(value.unique_id(), WORDS[1]);
    let mut plain = 222u64;
    value.set_unique_id(111);
    value.swap_unique_id(&mut plain);
    assert_eq!(value.unique_id(), 222);
    assert_eq!(plain, 111);
    *value.unique_id_mut() = 0;
    assert!(value.has_invalid_unique_id());
}

#[test]
fn source_text_cases_and_defined_unsigned_overflow() {
    let mut value = UniqueId(1_000_000);
    for (text, expected) in [
        ("", 0),
        ("_", 0),
        ("17", 17),
        ("18", 18),
        ("asdf_19", 19),
        ("_20", 20),
        ("_021", 21),
        ("asdf_19_22", 22),
        ("_20_23", 23),
        ("_021_024", 24),
        ("   _021_025     ", 0),
        ("bla", 0),
        ("123 456", 0),
        ("123_ 456", 0),
        ("123 _456", 456),
        ("123 _ 456", 0),
        ("123 456_", 0),
        ("123 456_  ", 0),
        ("123 456 _ ", 0),
        ("123 456  _", 0),
        ("123 456 ff", 0),
        ("_021bla_", 0),
        ("_021 bla    ", 0),
        ("_021 bla", 0),
        ("_021_bla", 0),
        ("_021_bla_", 0),
        ("18446744073709551615", u64::MAX),
        ("18446744073709551616", 0),
        ("18446744073709551617", 1),
        ("-1", 0),
        ("+1", 0),
        ("１２", 0),
        ("é_42", 42),
    ] {
        value.set_unique_id(1_000_000);
        value.set_unique_id_from_str(text).unwrap();
        assert_eq!(value.unique_id(), expected, "{text:?}");
    }
    value.set_unique_id(17);
    assert!(
        value
            .set_unique_id_from_str(&"0".repeat(1024 * 1024 + 1))
            .is_err()
    );
    assert_eq!(value.unique_id(), 17);
    value
        .set_unique_id_from_str(&"0".repeat(1024 * 1024))
        .unwrap();
    assert_eq!(value.unique_id(), 0);
}

#[test]
fn source_large_seed_replay_and_distinct_sample() {
    let mut random = UniqueIdGenerator::from_seed(SEED);
    let words: Vec<_> = (0..100_000).map(|_| random.get_unique_id()).collect();
    assert_eq!(
        words.iter().copied().collect::<HashSet<_>>().len(),
        words.len()
    );
    random.set_seed(SEED);
    for expected in words {
        assert_eq!(random.get_unique_id(), expected);
    }
}

#[test]
fn caller_owned_mutex_sharing_matches_one_serial_stream() {
    use std::sync::{Arc, Mutex};
    let shared = Arc::new(Mutex::new(UniqueIdGenerator::from_seed(SEED)));
    let threads: Vec<_> = (0..4)
        .map(|_| {
            let shared = Arc::clone(&shared);
            std::thread::spawn(move || {
                (0..250)
                    .map(|_| shared.lock().unwrap().get_unique_id())
                    .collect::<Vec<_>>()
            })
        })
        .collect();
    let mut actual: Vec<_> = threads
        .into_iter()
        .flat_map(|thread| thread.join().unwrap())
        .collect();
    let mut random = UniqueIdGenerator::from_seed(SEED);
    let mut expected: Vec<_> = (0..1000).map(|_| random.get_unique_id()).collect();
    actual.sort_unstable();
    expected.sort_unstable();
    assert_eq!(actual, expected);
}
