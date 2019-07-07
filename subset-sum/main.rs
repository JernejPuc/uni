mod subset_sum;
use std::{
    env,
    fs::{File, OpenOptions},
    io::{BufRead, BufReader, Write},
    str::FromStr,
    time::{Duration, Instant}
};


fn usize_from_line(buffer: &mut BufReader<File>) -> usize {
    // Read line of data (single int), convert to usize
    let mut tmp = String::new();
    let _ = buffer.read_line(&mut tmp);
    let _ = tmp.pop();
    let _ = tmp.pop();
    
    usize::from_str(&tmp).unwrap()
}


fn time_perf(f: &Fn(&Vec<usize>, usize, usize, f64) -> usize,
             a: &Vec<usize>, n: usize, k: usize, eps: f64) -> (f64, usize) {
    // Time performance of a chosen algorithm
    let now = Instant::now();
    let result: usize = f(&a, n, k, eps);
    let elapsed = now.elapsed();
    let dt: f64 =
        elapsed.as_secs() as f64
        + elapsed.subsec_millis() as f64 / 1000.0
        + elapsed.subsec_micros() as f64 / 1000000.0;

    (dt, result)
}


fn main() -> std::io::Result<()> {
    // CLI args
    let args: Vec<String> = env::args().collect();
    let alg: &str = args[1].trim();
    let infile: &str = args[2].trim();
    let mut eps: f64 = 0.5;

    if args.len() == 5 {
        eps = f64::from_str(&args[4]).unwrap();
    } else if args.len() == 4 && alg == "fptas" {
        eps = f64::from_str(&args[3]).unwrap();   
    }

    // Read set data
    let n: usize;
    let k: usize;
    let mut a: Vec<usize>;

    {
        let instance = File::open(infile)?;
        let mut buffer = BufReader::new(instance);

        n = usize_from_line(&mut buffer);
        k = usize_from_line(&mut buffer);
        a = vec![Default::default(); n];

        for i in 0..n {
            a[i] = usize_from_line(&mut buffer);
        }
    }

    // Execute
    let dt: f64;
    let result: usize;
    let multires: (f64, usize);

    match alg {
        "dynamic" => multires = time_perf(&subset_sum::dynamic, &a, n, k, 0.0),
        "exhaustive" => multires = time_perf(&subset_sum::exhaustive, &a, n, k, 0.0),
        "greedy" => multires = time_perf(&subset_sum::greedy, &a, n, k, 0.0),
        "fptas" => multires = time_perf(&subset_sum::fptas, &a, n, k, eps),
        _ => panic!("Algorithm '{}' not implemented.", alg)
    }

    dt = multires.0;
    result = multires.1;
    
    println!("{}:\nn = {}\nk = {}\nr = {}\nt = {}", alg, n, k, result, dt);
    
    // Write results
    if args.len() == 5 || (args.len() == 4 && alg != "fptas") {
        let outfile: &str = args[3].trim();

        let mut file = OpenOptions::new()
        .append(true)
        .create(true)
        .open(outfile)
        .unwrap();

        file.write((format!("{}:\nn = {}\nk = {}\nr = {}\nt = {}\n",
                            alg, n, k, result, dt)).as_bytes())?;

        if alg == "fptas" {
            file.write((format!("e = {}\n", eps)).as_bytes())?;
        }

        file.write(b"\n")?;
    }

    // End
    Ok(())
}
