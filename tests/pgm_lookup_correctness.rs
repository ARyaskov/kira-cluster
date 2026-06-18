#[test]
fn kira_pgm_finds_present_keys() {
    let keys: Vec<u64> = (0..10_000u64).map(|v| v * 3 + 7).collect();
    let pgm = kira_kv_engine::PgmIndex::build(keys.clone(), 64).unwrap();

    for (idx, key) in keys.iter().enumerate().step_by(113) {
        assert_eq!(pgm.index(*key).unwrap(), idx);
    }

    let absent = 8u64;
    assert!(pgm.index(absent).is_err());
}
