use enum_variants_strings::EnumVariantsStrings;
use ::std::ffi::c_int;


#[derive(Clone, Copy, EnumVariantsStrings)]
#[enum_variants_strings_transform(transform="none")]
pub enum Medium {
    InnerCore,
    OuterCore,
    LowerMantle,
    Mantle2,
    Mantle1,
    Mantle0,
    UpperMantle,
    LowerCrust,
    UpperCrust,
    Ocean,
    Troposphere0,
    Troposphere1,
    Stratosphere,
    Mesosphere,
    Exosphere,
    Exit,
    Topography,
    Unknown,
}

impl From<(c_int, bool)> for Medium {
    fn from(value: (c_int, bool)) -> Self {
        let (medium, ocean) = value;
        match medium {
            0 => Self::InnerCore,
            1 => Self::OuterCore,
            2 => Self::LowerMantle,
            3 => Self::Mantle2,
            4 => Self::Mantle1,
            5 => Self::Mantle0,
            6 => Self::UpperMantle,
            7 => Self::LowerCrust,
            8 => Self::UpperCrust,
            9 => if ocean { Self::Ocean } else { Self::UpperCrust },
            10 => Self::Troposphere0,
            11 => Self::Troposphere1,
            12 => Self::Stratosphere,
            13 => Self::Mesosphere,
            14 => Self::Exosphere,
            -1 => Self::Exit,
            100 => Self::Topography,
            _ => Self::Unknown,
        }
    }
}

impl From<Medium> for &'static str {
    fn from(value: Medium) -> Self {
        value.to_str()
    }
}

impl ::std::fmt::Display for Medium {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> ::std::fmt::Result {
        write!(f, "{}", Into::<&str>::into(*self))
    }
}

impl Medium {
    pub fn is_atmosphere(&self) -> bool {
        match self {
            Self::Troposphere0 => true,
            Self::Troposphere1 => true,
            Self::Stratosphere => true,
            Self::Mesosphere => true,
            Self::Exosphere => true,
            _ => false,
        }
    }
}
