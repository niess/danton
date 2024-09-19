use crate::bindings::danton;
use pyo3::prelude::*;
use pyo3::create_exception;
use pyo3::exceptions::{
    PyException, PyFileNotFoundError, PyIndexError, PyIOError, PyKeyboardInterrupt, PyKeyError,
    PyMemoryError, PyNotImplementedError, PySystemError, PyTypeError, PyValueError
};
use pyo3::ffi::PyErr_CheckSignals;
use ::std::ffi::{c_int, CStr};
use ::std::ptr::null_mut;


// ===============================================================================================
//
// Normalised errors.
//
// ===============================================================================================

#[derive(Debug, Default)]
pub struct Error<'a> {
    pub kind: Option<ErrorKind>,
    pub what: Option<&'a str>,
    pub why: Option<&'a str>,
}

#[derive(Clone, Copy, Debug)]
pub enum ErrorKind {
    CLibraryException,
    Exception,
    FileNotFoundError,
    IndexError,
    IOError,
    KeyboardInterrupt,
    KeyError,
    MemoryError,
    NotImplementedError,
    SystemError,
    TypeError,
    ValueError,
}

impl<'a> Error<'a> {
    pub fn maybe_what(mut self, what: Option<&'a str>) -> Self {
        self.what = what;
        self
    }

    pub fn maybe_why(mut self, why: Option<&'a str>) -> Self {
        self.why = why;
        self
    }

    pub fn new(kind: ErrorKind) -> Self {
        Self {
            kind: Some(kind),
            what: None,
            why: None,
        }
    }

    pub fn to_err(&self) -> PyErr {
        self.into()
    }

    pub fn to_string(&self) -> String {
        self.into()
    }

    pub fn what(mut self, what: &'a str) -> Self {
        self.what = Some(what);
        self
    }

    pub fn why(mut self, why: &'a str) -> Self {
        self.why = Some(why);
        self
    }
}

impl<'a> From<&Error<'a>> for String {
    fn from(value: &Error<'a>) -> Self {
        let Error { what, why, .. } = value;
        match what {
            None => match why {
                None => "something bad happened".to_string(),
                Some(why) => why.to_string(),
            }
            Some(what) => match why {
                None => format!("bad {what}"),
                Some(why) => format!("bad {what} ({why})"),
            },
        }
    }
}

impl<'a> From<Error<'a>> for String {
    fn from(value: Error<'a>) -> Self {
        (&value).into()
    }
}

impl<'a> From<&Error<'a>> for PyErr {
    fn from(value: &Error<'a>) -> Self {
        let msg: String = value.into();
        let kind = value.kind
            .unwrap_or(ErrorKind::Exception);
        match kind {
            ErrorKind::CLibraryException => PyErr::new::<CLibraryException, _>(msg),
            ErrorKind::Exception => PyErr::new::<PyException, _>(msg),
            ErrorKind::FileNotFoundError => PyErr::new::<PyFileNotFoundError, _>(msg),
            ErrorKind::IndexError => PyErr::new::<PyIndexError, _>(msg),
            ErrorKind::IOError => PyErr::new::<PyIOError, _>(msg),
            ErrorKind::KeyboardInterrupt => PyErr::new::<PyKeyboardInterrupt, _>(msg),
            ErrorKind::KeyError => PyErr::new::<PyKeyError, _>(msg),
            ErrorKind::MemoryError => PyErr::new::<PyMemoryError, _>(msg),
            ErrorKind::NotImplementedError => PyErr::new::<PyNotImplementedError, _>(msg),
            ErrorKind::SystemError => PyErr::new::<PySystemError, _>(msg),
            ErrorKind::TypeError => PyErr::new::<PyTypeError, _>(msg),
            ErrorKind::ValueError => PyErr::new::<PyValueError, _>(msg),
        }
    }
}

impl<'a> From<Error<'a>> for PyErr {
    fn from(value: Error<'a>) -> Self {
        (&value).into()
    }
}


// ===============================================================================================
//
// Variants explainers.
//
// ===============================================================================================

pub fn variant_explain(value: &str, options: &[&str]) -> String {
    let n = options.len();
    let options = match n {
        0 => unimplemented!(),
        1 => format!("'{}'", options[0]),
        2 => format!("'{}' or '{}'", options[0], options[1]),
        _ => {
            let options: Vec<_> = options
                .iter()
                .map(|e| format!("'{}'", e))
                .collect();
            format!(
                "{} or {}",
                options[0..(n - 1)].join(", "),
                options[n - 1],
            )
        },
    };
    format!(
        "expected one of {}, found '{}'",
        options,
        value
    )
}


// ===============================================================================================
//
// Danton interface.
//
// ===============================================================================================

create_exception!(danton, CLibraryException, PyException, "A C-library exception.");

pub fn to_result(rc: c_int, context: Option<*mut danton::Context>) -> Result<(), PyErr> {
    if rc == danton::SUCCESS {
        Ok(())
    } else {
        let context = context.unwrap_or_else(|| null_mut());
        let why = unsafe {
            let nerr = danton::error_count(context);
            if nerr == 0 {
                "an unexpected error occured".to_string()
            } else {
                let why = CStr::from_ptr(danton::error_pop(context))
                    .to_str()?
                    .to_string();
                for _ in 1..nerr {
                    danton::error_pop(context);
                }
                why
            }
        };
        let err = Error::new(ErrorKind::CLibraryException).why(&why);
        Err(err.to_err())
    }
}


// ===============================================================================================
//
// Keyboard interupts (catched by Python runtime).
//
// ===============================================================================================

pub fn ctrlc_catched() -> bool {
    if unsafe { PyErr_CheckSignals() } == -1 { true } else {false}
}
