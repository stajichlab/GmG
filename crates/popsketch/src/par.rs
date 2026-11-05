//! Parallel iteration that degrades to a single thread: with the `parallel` feature this is rayon's prelude; without it
//! `par_iter` and `into_par_iter` are the ordinary iterators (the browser build has no threads).

#[cfg(feature = "parallel")]
pub use rayon::prelude::*;

#[cfg(not(feature = "parallel"))]
pub trait ParIter {
    type Item;
    fn par_iter(&self) -> std::slice::Iter<'_, Self::Item>;
}

#[cfg(not(feature = "parallel"))]
impl<T> ParIter for [T] {
    type Item = T;
    fn par_iter(&self) -> std::slice::Iter<'_, T> {
        self.iter()
    }
}

#[cfg(not(feature = "parallel"))]
impl<T> ParIter for Vec<T> {
    type Item = T;
    fn par_iter(&self) -> std::slice::Iter<'_, T> {
        self.iter()
    }
}

#[cfg(not(feature = "parallel"))]
pub trait IntoParIter: Sized {
    type Iter: Iterator;
    fn into_par_iter(self) -> Self::Iter;
}

#[cfg(not(feature = "parallel"))]
impl IntoParIter for std::ops::Range<usize> {
    type Iter = std::ops::Range<usize>;
    fn into_par_iter(self) -> Self::Iter {
        self
    }
}
