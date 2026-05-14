import org.nlogo.build.{ NetLogoExtension, ExtensionDocumentationPlugin }

enablePlugins(NetLogoExtension, ExtensionDocumentationPlugin)

name := "transferentropy"
version := "1.0"
isSnapshot := true

scalaVersion := "3.7.0"
Compile / scalaSource := baseDirectory.value / "src" / "main"
scalacOptions ++= Seq("-deprecation", "-unchecked", "-Xfatal-warnings", "-feature", "-encoding", "us-ascii", "-release", "11")

netLogoVersion      := "7.0.3"
netLogoClassManager := "org.nlogo.extensions.transferentropy.TransferEntropyExtension"

libraryDependencies ++= Seq(
  "com.typesafe"       % "config"      % "1.3.1" % "test"
)
