
// Usage: Changes a list into its equivalent array form and vice versa.
function list_to_array(L) {
    const A = [];
    let i = 0;
    for (let p = L; !is_null(p); p = tail(p)) {
        A[i] = head(p);
        i = i + 1;
    }
    return A;
}

function array_to_list(a) {
    return build_list(i => a[i], array_length(a));
}

// Make a complete clone of any list / array based data structure.
// So you can wreck havoc on the thing without mutating the original data.
// Example:
// const obj = [list(1, list(2, 3), "string", [8, 9, list(3, 4)]);]
// const new_obj = deep_clone(obj);

function deep_clone(x) {
    if (is_pair(x)) {
        return pair(deep_clone(head(x)), deep_clone(tail(x)));
    } else if (is_array(x)) {
        const a = [];
        for (let i = 0; i < array_length(x); i = i + 1) {
            a[i] = deep_clone(x[i]);
        }
        return a;
    } else {
        return x;
    }
}

// Functions used to make debugging easier.
// Usage:
// Equivalent to display().
// Adapts itself to display_list() automatically if the input is a list.

function dynamic_display(...args) {
    // usage check
    if (array_length(args) === 0 || array_length(args) > 2) {
        error("dynamic_display(): Only 1 or 2 arguments are allowed.");
    }
    
    const x = args[0];
    const disp_fn = is_list(x) ? display_list : display;
    
    if (array_length(args) === 2) {
        disp_fn(x, args[1]);
    } else {
        disp_fn(x);
    }
}

const dd = dynamic_display;

// Usage:
// dbg(var1, var2, var3, ...)

// Example:
// dbg(math_PI, math_E);
function dbg(...args) {
    for (let i = 0; i < array_length(args); i = i + 1) {
        dynamic_display(args[i]);
    }
}

// Usage:
// dbgx(()=>var1, ()=>var2, ()=>var3, ...)

// Example:
// dbgx(()=>math_PI, ()=>math_E);
function dbgx(...args) {
    // usage check
    for (let i = 0; i < array_length(args); i = i + 1) {
        if (!is_function(args[i]) || arity(args[i]) !== 0) {
            error("dbgx(): Arguments must be nullary functions");
        }
    }
    
    // helper fns to get var name
    const name_of = f => stringify(f);
    
    function substr(s, n) {
        let new_str = "";
        for (let i = n; i < array_length(s); i = i + 1) {
            new_str = new_str + char_at(s, i);
        }
        return new_str;
    }
    
    const var_name_of = f => substr(name_of(f), 6);
    
    // print
    for (let i = 0; i < array_length(args); i = i + 1) {
        const curr_arg = args[i];
        
        dynamic_display(curr_arg(), var_name_of(curr_arg) + ":");
    }
}

// Functions for simple data structures.
// Maintains a stack data structure, first in last out.
const stack = [];
let stack_length = 0;

function push(x) {
    stack[stack_length] = x;
    stack_length = stack_length + 1;
}

function pop() {
    // usage check
    if (stack_length === 0) {
        error("pop() cannot be done on a stack of size 0.");
    }
    
    stack_length = stack_length - 1;
    return stack[stack_length];
}

// Maintains a queue data structure, first in first out.
let queue = list();

function enqueue(x) {
    queue = append(queue, list(x));
}

function dequeue() {
    if (is_null(queue)) {
        error("dequeue() cannot be done on a queue of size 0.");
    }
    
    const ans = head(queue);
    queue = tail(queue);
    return ans;
}

// Semi-useful pre-supplied functions

// Permutations
// Input: List
// Output: List of all permutations
// https://github.com/source-academy/source-programs

// Usage:
// permutations(list(1, 2, 3));
function flatmap(f, seq) {
    return accumulate(append, null, map(f, seq));
}

function permutations(s) {
    return is_null(s)             // empty set?
           ? list(null)           // sequence containing empty set
           : flatmap(x => map(p => pair(x, p),
                              permutations(remove(x, s))),
                     s);
}


// Subsets
// Input: List
// Output: List of all subsets
// https://github.com/source-academy/source-programs

// Usage:
// subsets(list(1, 2, 3));
function subsets(s) {
    if (is_null(s)) {
        return list(null);
    } else {
        const rest = subsets(tail(s));
        return append(rest, map(x => pair(head(s), x), rest));
    }
}

// Merge Sort
// Input: List
// Output: Sorted List
// https://github.com/source-academy/source-programs
function merge_sort(xs) {
    if ( is_null(xs) || is_null(tail(xs))) {
        return xs;
    } else {
        const mid = middle(length(xs)) ;
        return merge(merge_sort(take(xs, mid)),
                     merge_sort(drop(xs, mid)));
    }
}

function middle(n) {
    return math_floor(n / 2) ;
}

function merge(xs , ys) {
    if ( is_null(xs)) {
        return ys;
    } else if (is_null(ys)) {
        return xs;
    } else {
        const x = head(xs) ;
        const y = head(ys) ;
        if (x <= y) {
            return pair(x, merge(tail(xs), ys));
        } else {
            return pair(y, merge(xs, tail(ys)));
        }
    }
}

function take(xs , n) {
    return (n === 0) 
        ? null
        : pair(head(xs), take(tail(xs), n - 1));
}

function drop(xs, n) {
    return n === 0
        ? xs
        : drop(tail(xs), n - 1);
}

// Bubble Sort
// Input: Array
// Output: Sorted Array * - sorts in place
// https://www.crio.do/blog/top-10-sorting-algorithms/
function array_sort(a) {
    function swap(a, i, j) {
        const tmp = a[i];
        a[i] = a[j];
        a[j] = tmp;
    }
    
    const n = array_length(a);
    for (let i = 0; i < n - 1; i = i + 1) {
        for (let j = 0; j < n - i - 1; j = j + 1) {
            if (a[j] > a[j + 1]) {
                swap(a, j, j + 1);
            }
        }
    }
    
    return a;
}



function rotate_matrix(M){
    const l= array_length(M);
    const ll= math_floor(l/2);
    
    for (let i=0; i<l; i=i+1){
        for (let j=0; j<i; j=j+1){
            let temp= M[i][j];
            M[i][j]= M[j][i];
            M[j][i]=temp;
        }
    }
    
    for (let i=0; i<l; i=i+1){
        for (let j=0; j<ll; j=j+1){
            let temp= M[i][j];
            M[i][j]= M[i][l-j-1];
            M[i][l-j-1]=temp;
        }
    }
    
}
let A = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16]];
rotate_matrix(A);
A;


function stream_pairs(s) {
return is_null(s)
? null
: stream_append(
stream_map(
sn => pair(head(s), sn),
stream_tail(s)),
stream_pairs(stream_tail(s)));
}

ints = streams(1, 2, 3, 4, 5);

function accumulate_tree(f, op, initial, tree) {
return accumulate(
(x, ys) => !is_list(x)
? op(f(x), ys)
: op(accumulate_tree(f, op, initial, x), ys),
initial,
tree );
}

function count_data_items(tree) {
return is_null(tree)
? 0
: ( is_list(head(tree))
? count_data_items(head(tree))
: 1 )
+
count_data_items(tail(tree));

function scale_tree(tree, k) {
return map(sub_tree =>
!is_list(sub_tree)
? k * sub_tree
: scale_tree(sub_tree, k),
tree);

function partial_sums(s) {
    let acc = 0;
    // YOUR SOLUTION HERE
    function helper(s) {
        acc = head(s) + acc;
        return pair(acc, 
                    () => helper(stream_tail(s)));
    }
    return helper(s);
}

function partial_sums2(s) {
    return pair(head(s), 
    () => stream_map(x => x + head(s), partial_sums2(stream_tail(s))));

}

function sum_odd(n) {
    const term = x => x;
    const next = x => x + 2;
    return sum(term, 1, next, n * 2);


function sum(term, a, next, b) {
    return a > b
      ? 0
      : term(a) + sum(term, next(a), next, b);
}

function prime_only(xs) {
    function check(n) {
        
        function check_helper(n, d) {
            return n === d
            ? true
            : n % d === 0
            ? false
            : check_helper(n, d + 1);
        }
    
    return n === 1 ? false : check_helper(n, 2);
    
    }  
    
    return filter( check , xs);
}



function shorten_stream(s, k) {
    if (k === 0) {
        return null;
    }
    else if ( is_null(s)) { // k >= length of s
        return s;
    }
    else {
        return pair(head(s), () => shorten_stream(stream_tail(s), k - 1));
    }
}


function alternating() {
  return pair(1,
                () => pair(-1, 
                            () => alternating()));
}

const alternating_ones = alternating();


function make_alternating_stream(s) {
    return pair(head(s),
            () => pair(-head(stream_tail(s)),
            () => make_alternating_stream(stream_tail(stream_tail(s)))));
}

function zip_streams(s1, s2) {
    return is_null(s1)
    ? s2
    :is_null(s2)
    ? s1
    : pair(head(s1), () => pair(head(s2), 
                     () => zip_streams(stream_tail(s1), stream_tail(s2))));
}

function every_other(s) {
    return pair(head(s), () => every_other(stream_tail(stream_tail(s))));
}


function partial_sums(s) {
    return pair(head(s), 
    () => stream_map(x => x + head(s), partial_sums(stream_tail(s))));

}

// TASK 1

function d_split_list(xs) {
    let pointer = xs;
  const middle = math_ceil(length(xs)/2);
  for (let i = 1 ; i < middle ; i = i + 1) {
      pointer = tail(pointer);
  }
  const back = tail(pointer);
  set_tail(pointer, null);
  return pair(xs, back);
}

// TEST:
// const my_list1 = list(1, 2, 3, 4, 5, 6);
// const my_list2 = list(5, 4, 3, 2, 1);
// d_split_list(my_list1);
// d_split_list(my_list2);

// TASK 2
 

function d_merge(xs, ys) {
    if (is_null(xs)) {
        return ys;
    }
    else if (is_null(ys)) {
        return xs;
    }
    else {
        if (head(xs) < head(ys)) {
            set_tail(xs, d_merge(tail(xs), ys));
            return xs;
            
        } else {
            set_tail(ys, d_merge(xs, tail(ys)));
            return ys;
        }
    }
}

// TEST:
const my_list1 = list(2, 4, 5, 9);
const my_list2 = list(3, 5, 8);
d_merge(my_list1, my_list2);


function d_merge_sort(xs) {
    if (is_null(xs) || is_null(tail(xs))) {
        return xs;
    } else {
        const split = d_split_list(xs);
        return d_merge(d_merge_sort(head(split)), d_merge_sort(tail(split)));
    }

}

// TEST:
const my_list = list(7, 2, 4, 6, 9, 1, 5, 8, 3, 6);
d_merge_sort(my_list);


function is_entry(key, dict) {
return accumulate( (x, y) => (head(x) === key) || y
,
false ,
dict );
}

function remove_entry(key, dict) {
return filter( x => head(x) !== key
,
dict );
}

function circular_left_shift(xs) {
if (is_null(xs)) {
return xs;
} else {
// Solution 1:
return reverse(pair(head(xs), reverse(tail(xs))));
// Solution 2:
return accumulate(pair, list(head(xs)), tail(xs));
// Solution 3:
return accumulate((x, y) => is_null(y)
? list(x, head(xs)) : pair(x, y),
null, tail(xs));
// Solution 4:
return is_null(tail(xs))
? xs
: pair(head(tail(xs),
circular_left_shift(pair(head(xs),
tail(tail(xs))))));
}
}

function circular_right_shift(xs) {
if (is_null(xs)) {
return xs;
} else {
// Solution 1:
const rev = reverse(xs);
return pair(head(rev), reverse(tail(rev)));
// Solution 2:
if (is_null(tail(xs))) {
return xs;
} else {
const wish = circular_right_shift(tail(xs));
return pair(head(wish), pair(head(xs), tail(wish)));
}
}
}

accumulate((x, ys) => circular_left_shift(pair(x, ys)), null, xs);

function max_tree(tree) {
return is_null(tree)
? 0
: is_list(head(tree))
? math_max(max_tree(head(tree)), max_tree(tail(tree)))
: math_max(head(tree), max_tree(tail(tree)));
}

function alt_column_matrix(R, C) {
    // Initialise the whole r by c matrix
    const res = [];
    for (let i = 0; i < R; i = i + 1) {
        res[i] = [];
    }
    
    function next(r, c) {
        if (c % 2 === 0 && r < R - 1) {
            return pair(r + 1, c);
        }
        if (c % 2 === 0 && r === R - 1) {
            return pair(r, c + 1);
        }
        if (c % 2 === 1 && r > 0) {
            return pair(r - 1, c);
        }
        if (c % 2 === 1 && r === 0) {
            return pair(r, c + 1);
        }
    }
    
    const total_iterations = R * C;
    let curr = pair(0, 0);
    for (let i = 1; i <= total_iterations; i = i + 1) {
        const r = head(curr);
        const c = tail(curr);
        res[r][c] = i;
        curr = next(r, c);
    }
    return res;
}












// ========================================================
// MY OWN STUFF 
// ========================================================
// LIST TO ARRAY
// ========================================================
function list_to_array(xs){
    const res=[];
    let current_pair= xs;
    while(!is_null(current_pair)){
        res[array_length(res)]= head(current_pair);
        current_pair = tail(current_pair);
    }
    return res;
}

// ========================================================
// ARRAY TO LIST
// ========================================================
function array_to_list(a){
    let res=null;
    const len= array_length(a);
    for(let i=1; i<= len; i=i+1){
        res= pair(a[len-i], res);
    }
    return res;
}

// ========================================================
// MERGE SORT (LISTS)
// ========================================================
function middle(n) {
    return math_floor(n / 2);
}

// put the first n elements of xs into a list
function take(xs, n) {
    return n === 0
    ? null 
    : pair(head(xs), take(tail(xs), n - 1));
}

// drop the first n elements from list, return rest
function drop(xs, n) {
    return n === 0
    ? xs
    : drop(tail(xs), n - 1);
}

function merge(xs, ys) {
    if (is_null(xs)) {
        return ys;
    } else if (is_null(ys)) {
        return xs;
    } else {
        const x = head(xs);
        const y = head(ys);
        return x < y
            ? pair(x, merge(tail(xs), ys))
            : pair(y, merge(xs, tail(ys)));
    }
}

function merge_sort(xs) {
    if (is_null(xs) || is_null(tail(xs))) {
        return xs;
    } else {
        const mid = middle(length(xs));
        return merge(merge_sort(take(xs, mid)),
        merge_sort(drop(xs, mid)));
    }
}

merge_sort(list(10, 8, 2, 5, 6, 1));

// ========================================================
// MERGE SORT (ARRAYS)
// ========================================================
function merge_sort(A) {
    merge_sort_helper(A, 0, array_length(A) - 1);
}

function merge_sort_helper(A, low, high) {
    if (low < high) {
        const mid = math_floor((low + high) / 2);
        merge_sort_helper(A, low, mid);
        merge_sort_helper(A, mid + 1, high);
        merge(A, low, mid, high);
    }
}

function merge(A, low, mid, high) {
    const B = [];
    let left = low;
    let right = mid + 1;
    let Bidx = 0;
    
    while (left <= mid && right <= high) {
        if (A[left] <= A[right]) {
            B[Bidx] = A[left];
            left = left + 1;
        } else {
            B[Bidx] = A[right];
            right = right + 1;
        }
        Bidx = Bidx + 1;
    }
    
    while (left <= mid) {
        B[Bidx] = A[left];
        Bidx = Bidx + 1;
        left = left + 1;
    }   
    while (right <= high) {
        B[Bidx] = A[right];
        Bidx = Bidx + 1;
        right = right + 1;
    }
    
    for (let k = 0; k < high - low + 1; k = k + 1) {
        A[low + k] = B[k];
    }
}

const A = [3, 9, 2, 1, 6, 5, 3, 8];
merge_sort(A);
A;

// ========================================================
// QUICK SORT (LISTS)
// ========================================================

function partition(xs, p) {
    function helper(lst, less, more) {
        return is_null(lst)
        ? pair(less, more)
        : head(lst) <= p
        ? helper(tail(lst), pair(head(lst), less), more)
        : helper(tail(lst), less, pair(head(lst), more)); //head(lst) > p
        }
        return helper(xs, null, null);
}

function quicksort(xs) {
    if (is_null(xs) || is_null(tail(xs))) {
        return xs;
    } else {
    const pivot = head(xs);
    const smaller = head(partition(tail(xs), head(xs)));
    const bigger = tail(partition(xs, head(xs)));
    return append(quicksort(smaller), pair(pivot, quicksort(bigger)));
    }
}

// Test
const my_test_list = list(23, 12, 56, 92, -2, 0);
quicksort(my_test_list); 

// ========================================================
// INSERTION SORT (LISTS)
// ========================================================

function insert(x, xs) {
    return is_null(xs)
        ? list(x)
        : x <= head(xs)
        ? pair(x, xs)
        : pair(head(xs), insert(x, tail(xs)));
}

function insertion_sort(xs) {
    return is_null(xs)
        ? xs
        : insert(head(xs), insertion_sort(tail(xs)));
}

const my_test_list = list(23, 12, 56, 92, -2, 0);
insertion_sort(my_test_list);

// ========================================================
// INSERTION SORT (ARRAYS)
// ========================================================

function insertion_sort(A) {
    const len = array_length(A);
    
    for (let i = 1; i < len; i = i + 1) {
        let j = i - 1;
        while (j >= 0 && A[j] > A[j + 1]) {
            swap(A, j, j + 1);
            j = j - 1;
        }
    }
}

function swap(A, x, y) {
    const temp = A[x];
    A[x] = A[y];
    A[y] = temp;
}

const A = [3, 9, 2, 1, 6, 5, 3, 8];
insertion_sort(A);
A;

// ========================================================
// GENERALISED INSERTION SORT (LISTS)
// ========================================================

function insert_cmp(x, xs, cmp) {
    return is_null(xs) 
          ? list(x)
          : cmp(x, head(xs)) 
          ? pair(x, xs)
          : pair(head(xs), insert_cmp(x, tail(xs), cmp));
}

function insertion_sort_cmp(xs, cmp) {
    return is_null(xs) 
          ? xs
          : insert_cmp(head(xs), 
                        insertion_sort_cmp(tail(xs), cmp),
                        cmp);
}

// Test
const xs = list(6, 3, 8, 5, 1, 9, 6, 4, 2, 7);

function is_even(x) {
    return x % 2 === 0;
}

function is_odd(x) {
    return !is_even(x);
}

(a)
insertion_sort_cmp(xs, (x, y) => y > x);
Result: list(1, 2, 3, 4, 5, 6, 6, 7, 8, 9)

(b)
insertion_sort_cmp(xs, (x, y) => x > y);
Result: list(9, 8, 7, 6, 6, 5, 4, 3, 2, 1)

(c)
insertion_sort_cmp(xs, (x, y) => false);
Result: list(7, 2, 4, 6, 9, 1, 5, 8, 3, 6)

(d)
insertion_sort_cmp(xs, (x, y) => is_even(x) && is_even(y)
                                ? x < y
                                : is_even(x) && is_odd(y)
                                ? true 
                                : is_odd(x) && is_even(y)
                                ? false
                                : is_odd(x) && is_odd(y)
                                ? x > y
                                : false);
Result: list(2, 4, 6, 6, 8, 9, 7, 5, 3, 1)

// ========================================================
// SELECTION SORT (LISTS)
// ========================================================

function smallest(xs) {
    return accumulate((x, y) => x < y ? x : y,
        head(xs), tail(xs));
}

function selection_sort(xs) {
    if (is_null(xs)) {
        return xs;
    } else {
        const x = smallest(xs);
        return pair(x, selection_sort(remove(x, xs)));
    }
}

const my_test_list = list(23, 12, 56, 92, -2, 0);
selection_sort(my_test_list);

// ========================================================
// SELECTION SORT (ARRAYS)
// ========================================================

function selection_sort(A) {
    const len = array_length(A);

    for (let i = 0; i < len - 1; i = i + 1) {
        let min_pos = find_min_pos(A, i, len - 1);
        swap(A, i, min_pos);
    }
}

function find_min_pos(A, low, high) {
    let min_pos = low;
    for (let j = low + 1; j <= high; j = j + 1) {
        if (A[j] < A[min_pos]) {
            min_pos = j;
        }
    }
    return min_pos;
}

function swap(A, x, y) {
    const temp = A[x];
    A[x] = A[y];
    A[y] = temp;
}

const A = [3, 9, 2, 1, 6, 5, 3, 8];
selection_sort(A);
A;

// ========================================================
// BINARY SEARCH (ARRAYS)
// ========================================================

function binary_search(A, v) {
    let low = 0;
    let high = array_length(A) - 1;

    while (low <= high) {
        const mid = math_floor((low + high) / 2 );
        if (v === A[mid]) {
            break;
        } else if (v < A[mid]) {
            high = mid - 1;
        } else {
            low = mid + 1;
        }
    }
    return (low <= high);
}

binary_search([1,2,3,4,5,6,7,8,9], 8);

// ========================================================
// ARRAY PROCESSING
// ========================================================

function accumulate_array(op, init, A) {
    const len = array_length(A);
    let result = init;
    for (let i = 0; i < len; i = i + 1) {
        result = op(result, A[i]);
    }
    return result;
}

function filter_array(pred, A) {
    let result = [];
    let len = array_length(A);
    let j = 0;
    for (let i = 0 ; i < len ; i = i + 1) {
        if (pred(A[i])) {
            result[j] = A[i];
            j = j + 1;
        }
        
    }
    return result;
}

// ========================================================
// MEMOIZATION (READ WRITE IS HEREEEEEE)
// ========================================================

const mem = [];

function read(n, k) {
    return mem[n] === undefined
          ? undefined
          : mem[n][k];
}

function write(n, k, value) {
    if (mem[n] === undefined) {
        mem[n] = [];
    }
    mem[n][k] = value;
}

function mchoose(n, k) {
    if (read(n, k) !== undefined) {
        return read(n, k);
    } else {
        const result = k > n
                      ? 0
                      : k === 0 || k === n
                      ? 1
                      : mchoose(n - 1, k) + mchoose(n - 1, k - 1);
        write(n, k, result);
        return result;
    }
}

mchoose(12, 6);
mchoose(100, 50);

// ========================================================
// ACCUMULATE TREE
// ========================================================

function accumulate_tree(f, op, initial, tree) {
    return accumulate((curr, wish) => is_list(curr)
                ? op(accumulate_tree(f, op, initial, curr), wish)
                : op(f(curr), wish),
                initial,
                tree);
}

function tree_sum(tree) {
    return accumulate_tree(x => x, (x, y) => x + y, 0 , tree);
}

function count_data_items(tree) {
    return accumulate_tree(x => 1, (x, y) => x + y, 0 , tree);
}

function flatten(tree) {
    return accumulate_tree(x => list(x), append, null , tree);
}

// Test
const tree1 = list(1, 2, list(3, 4), list(5, list(6, 7)));

const tree2 = list(1, list(list(8, 9), 10, list(11, list(12))), 
                  null, list(3, 4), list(5, list(6, 7)));

display( tree_sum(tree1) ); // Result: 28
display( tree_sum(tree2) ); // Result: 76

display( count_data_items(tree1) ); // Result: 7
display( count_data_items(tree2) ); // Result: 11

display_list( flatten(tree1) );
// Result: list(1, 2, 3, 4, 5, 6, 7)
display_list( flatten(tree2) );
// Result: list(1, 8, 9, 10, 11, 12, 3, 4, 5, 6, 7)
// ========================================================
// CHECK CYCLE IN LINKED LIST
// ========================================================
function check_cycle(xs) {
    let slow = xs;
    let fast = xs;
    
    while (!is_null(fast) && !is_null(tail(fast))) {
        slow = tail(slow);
        fast = tail(tail(fast));
        
        if (slow === fast) {
            return true;
        }
    }
    
    return false;
}


function max(...args) {
    // Check if a single array is passed as an argument
    const values = array_length(args) === 1 && is_array(args[0]) ? args[0] : args;
    
    // Set the initial maximum to the first element
    let maxValue = values[0];
    
    // Loop through the array to find the maximum value
    for (let i = 1; i < array_length(values); i = i + 1) {
        if (values[i] > maxValue) {
            maxValue = values[i];
        }
    }
    
    return maxValue;
}


function min(...args) {
    // Check if a single array is passed as an argument
    const values = array_length(args) === 1 && is_array(args[0]) ? args[0] : args;
    
    // Set the initial minimum to the first element
    let minValue = values[0];
    
    // Loop through the array to find the minimum value
    for (let i = 1; i < array_length(values); i = i + 1) {
        if (values[i] < minValue) {
            minValue = values[i];
        }
    }
    
    return minValue;
}



Count Islands (BFS)
function count_islands(map) {
    let islands = 0;
    let visited = [];
    
    const rows = array_length(map);
    const cols = array_length(map[0]);
    
    for (let i = 0; i < rows; i = i+1) {
        visited[i] = [];
        for (let j = 0; j < cols; j = j+1) {
            visited[i][j] = 0;
        }
    }
    
    const dirs = [[1, 0], [0, 1], [-1, 0], [0, -1]];
    const checkBoundaryRow = r => r >= 0 && r < array_length(map);
    const checkBoundaryCol = c => c >= 0 && c < array_length(map[0]);

    
    function bfs(r,c) {
        let queue = list();

        function enqueue(x) {
            queue = append(queue, list(x));
        }
        
        function dequeue() {
            if (is_null(queue)) {
                error("dequeue() cannot be done on a queue of size 0.");
            }
    
            const ans = head(queue);
            queue = tail(queue);
            return ans;
        }
        
        visited[r][c] = 1;
        
        enqueue([r,c]);
        
        while (length(queue) > 0) {
            let point = dequeue();
            for (let i = 0; i < array_length(dirs); i = i+1) {
                let newR = point[0] + dirs[i][0];
                let newC = point[1] + dirs[i][1];
                if (checkBoundaryRow(newR) && checkBoundaryCol(newC) && map[newR][newC] !== 0 && visited[newR][newC] !== 1) {
                    enqueue([newR, newC]);
                    visited[newR][newC] = 1;
                }
            }
        }
    }
    for (let r = 0; r < array_length(map); r= r+1) {
        for (let c = 0; c < array_length(map[0]); c = c+1) {
            if (map[r][c] !== 0 && visited[r][c] !== 1) {
                bfs(r,c);
                islands = islands + 1;
            }
            
        }
    }
    
    return islands;
    
}
