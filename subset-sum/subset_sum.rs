pub fn dynamic(a: &Vec<usize>, n: usize, k: usize, _: f64) -> usize {
    let mut S: Vec<Vec<usize>> = vec![vec![0; k+1]; n+1];

    for i in 1..n+1 {
        for j in 1..k+1 {
            if a[i-1] > j {
                S[i][j] = S[i-1][j];
            } else {
                S[i][j] = std::cmp::max(S[i-1][j], S[i-1][j-a[i-1]] + a[i-1]);
            }
        }
    }

    S[n][k]
}


pub fn exhaustive(a: &Vec<usize>, n: usize, k: usize, _: f64) -> usize {
    let mut Li: Vec<usize> = vec![0];

    for i in 0..n {
        Li = merge_sort_and_filter(&Li, a[i], k);
    }

    Li.into_iter().max().unwrap()
}


pub fn greedy(a_: &Vec<usize>, n: usize, k: usize, _: f64) -> usize {
    let mut G: usize = 0;
    let mut a: Vec<usize> = a_.to_owned();
    a.sort_unstable_by(|a, b| b.cmp(a));

    for i in 0..n {
        if a[i] <= k - G {
            G += a[i];
        }
    }

    G
}


pub fn fptas(a: &Vec<usize>, n: usize, k: usize, eps: f64) -> usize {
    let mut Li: Vec<usize> = vec![0];
    let one_plus_delta = 1.0 + eps / ((2*n) as f64);

    for i in 0..n {
        Li = trim(&merge_sort_and_filter(&Li, a[i], k), one_plus_delta);
    }

    Li.into_iter().max().unwrap()
}


fn merge_sort_and_filter(a: &Vec<usize>, e: usize, k: usize) -> Vec<usize> {
    let b: Vec<usize> = a.iter().map(|x| x + e).collect();
    let a_len = a.len();
    let mut merged: Vec<usize> = vec![0; 2*a_len];
    let mut a_idx: usize = 0;
    let mut b_idx: usize = 0;
    
    for i in 0..merged.len() {
        if a_idx < a_len {
            if a[a_idx] < b[b_idx] {
                merged[i] = a[a_idx];
                a_idx += 1;
            } else if a[a_idx] == b[b_idx] {
                merged[i] = a[a_idx];
                a_idx += 1;
                b_idx += 1;
            } else {
                merged[i] = b[b_idx];
                b_idx += 1;
            }
        } else if b_idx < a_len && b[b_idx] <= k {
            merged[i] = b[b_idx];
            b_idx += 1;
            a_idx = i+1;
        } else {
            a_idx = i;
            break;
        }
    }

    merged.truncate(a_idx);
    merged
}


fn trim(a: &Vec<usize>, one_plus_delta: f64) -> Vec<usize> {
    let a_len = a.len();
    let mut last = a[0];
    let mut trimmed = vec![last; a_len];
    let mut idx = 1;

    for i in 1..a_len {
        if (a[i] as f64) > (last as f64) * one_plus_delta {
            trimmed[idx] = a[i];
            last = a[i];
            idx += 1;
        }
    }
    
    trimmed.truncate(idx);
    trimmed
}
