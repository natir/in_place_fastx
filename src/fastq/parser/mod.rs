/* mod declaration */
pub mod sequential;
pub mod shared_state;

pub use sequential::*;
pub use shared_state::*;

#[cfg(test)]
mod tests {
    pub type BaseCount<T> = [T; 4];
    pub trait AbsBaseCount {
        fn new() -> Self;
    }
}
