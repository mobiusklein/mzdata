//! Elements of controlled vocabularies used to describe mass spectra and their components.
//!
//! This crate implements the "CV param" model used throughout the PSI-MS family of mass
//! spectrometry file formats (mzML, imzML, mzMLb, ...): every piece of metadata is either
//! a term drawn from a controlled vocabulary (like the PSI-MS ontology) or a free-text,
//! user-defined key/value pair. Both are represented uniformly by [`Param`].
//!
//! ## Core types
//!
//! - [`Param`] / [`ParamCow`]: a single CV or user parameter. [`Param`] owns its data;
//!   [`ParamCow`] borrows where possible (useful for `const` tables of well-known terms).
//! - [`Value`] / [`ValueRef`]: the value half of a [`Param`], which may be a string,
//!   integer, float, boolean, byte buffer, list, or empty. [`Value`] owns its data;
//!   [`ValueRef`] borrows where possible. Both implement [`ParamValue`] for reading and
//!   coercing the stored value.
//! - [`CURIE`]: a `namespace:accession` identifier (e.g. `MS:1000041`) that names a term
//!   within a [`ControlledVocabulary`].
//! - [`ControlledVocabulary`]: the vocabulary a term belongs to (`MS`, `UO`, ...), and a
//!   factory for building [`Param`]s within that namespace.
//! - [`Unit`]: a closed set of units of measure (from the Unit Ontology and PSI-MS) that
//!   a [`Param`]'s value may be expressed in.
//! - [`ParamDescribed`] / [`ParamDescribedRead`]: traits for types that carry a list of
//!   [`Param`]s (spectra, scans, instrument components, ...), with helpers to look a
//!   parameter up by name, [`CURIE`], or accession string.
//!
//! ## Building a parameter
//!
//! Construct a controlled-vocabulary parameter from a [`ControlledVocabulary`] namespace:
//!
//! ```rust
//! use mzdata_param::{ControlledVocabulary, ParamLike, ParamValue};
//!
//! let param = ControlledVocabulary::MS.param_val("MS:1000041", "charge state", 2i32);
//! assert_eq!(param.name(), "charge state");
//! assert_eq!(param.curie().unwrap().to_string(), "MS:1000041");
//! assert_eq!(param.value().to_i64().unwrap(), 2);
//! ```
//!
//! Or use [`ParamBuilder`] / [`Param::builder`] for more incremental construction,
//! including plain user-defined parameters that have no controlled vocabulary:
//!
//! ```rust
//! use mzdata_param::{Param, ParamValue, Unit, curie};
//!
//! let p = Param::builder()
//!     .name("scan start time")
//!     .curie(curie!(MS:1000016))
//!     .value(12.34)
//!     .unit(Unit::Minute)
//!     .build();
//! assert_eq!(p.to_f64().unwrap(), 12.34);
//! assert_eq!(p.unit, Unit::Minute);
//!
//! // A user-defined parameter has no CURIE at all.
//! let custom = Param::new_key_value("my custom field", "some value");
//! assert!(!custom.is_controlled());
//! ```
//!
//! ## Reading values generically
//!
//! [`ParamValue`] lets you read a [`Param`]'s value without caring whether it was parsed
//! as a string, integer, or float - coercions are attempted on demand:
//!
//! ```rust
//! use mzdata_param::{ControlledVocabulary, ParamValue};
//!
//! let p = ControlledVocabulary::MS.param_val("MS:1000827", "isolation window target m/z", "500.25");
//! assert_eq!(p.to_f64().unwrap(), 500.25);
//! assert_eq!(p.to_str(), "500.25");
//! ```
//!
//! ## Attaching parameters to types with [`ParamDescribed`]
//!
//! The [`ParamDescribed`] trait provides many methods that make it easier to operate
//! on instances of types which are described by a list of [`Param`]s.
//!
//! ### Simple implementation helpers
//!
//! Use [`impl_param_described!`] (for a plain `Vec<Param>` field) or
//! [`impl_param_described_deferred!`] (for an `Option<Vec<Param>>` field that is
//! lazily allocated on first write) to implement [`ParamDescribed`] for your type:
//!
//! ```rust
//! use mzdata_param::{impl_param_described, Param, ParamDescribed, ParamList, ParamValue};
//!
//! #[derive(Default)]
//! struct MyComponent {
//!     params: ParamList,
//! }
//!
//! impl_param_described!(MyComponent);
//!
//! let mut c = MyComponent::default();
//! c.add_param(Param::new_key_value("vendor", "Acme"));
//! assert_eq!(c.get_param_by_name("vendor").unwrap().to_str(), "Acme");
//! ```
//!
//! ## Matching against known terms with `CURIE`
//!
//! [`CURIE`] implements `PartialEq` against anything implementing [`ParamLike`], so you can
//! compare a parameter directly against a well-known accession, and the [`curie!`] macro
//! gives a compact way to write one inline:
//!
//! ```rust
//! use mzdata_param::{ControlledVocabulary, ParamLike, curie};
//!
//! let p = ControlledVocabulary::MS.param("MS:1000016", "scan start time");
//! assert!(curie!(MS:1000016) == p);
//! ```
//!
//! ## The PSI-MS ontology (`cv` feature)
//!
//! With the `cv` feature enabled, [`MSVocabulary`] gives access to the full PSI-MS
//! ontology (term names, synonyms, and the `is_a` parent/child hierarchy), backed by an
//! embedded static snapshot and an optional on-disk cache that can be refreshed from a
//! `.obo` file via [`MSVocabulary::update_from_obo`]:
//!
//! ```rust
//! # #[cfg(feature = "static_data")]
//! # {
//! use mzdata_param::{MSVocabulary, curie};
//!
//! let term = MSVocabulary::get(curie!(MS:1000044)).unwrap();
//! assert_eq!(term.name.as_ref(), "dissociation method");
//! assert!(MSVocabulary::is_child_of(curie!(MS:1000133), curie!(MS:1000044)));
//! # }
//! ```
pub use mzdata_param::*;