package org.nlogo.extensions.transferentropy

import java.io

import org.nlogo.nvm.ExtensionContext

import scala.jdk.CollectionConverters.IteratorHasAsScala
import scala.reflect.Selectable.reflectiveSelectable
import org.nlogo.api.ScalaConversions.RichAny

import org.nlogo.core.LogoList
import org.nlogo.api._
import org.nlogo.core.{ NumberParser, Syntax }
import org.nlogo.api.Reporter
import org.nlogo.api.Command
import org.nlogo.core.Syntax._

import infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov
import infodynamics.utils.MatrixUtils


class TransferEntropyExtension extends DefaultClassManager {

  class TransferEntropy extends Reporter {
    override def getSyntax = Syntax.reporterSyntax(right = List(ListType, ListType, NumberType, NumberType), ret = NumberType)
    def report(args: Array[Argument], context: Context): AnyRef = {
      val X: Array[Double] = try args(0).getList.toArray.map(_.asInstanceOf[Double])
      catch {case e: LogoException => throw new ExtensionException(e.getMessage)}
      val Y: Array[Double] = try args(1).getList.toArray.map(_.asInstanceOf[Double])
      catch {case e: LogoException => throw new ExtensionException(e.getMessage)}

      val k = args(2).getIntValue
      val tau = args(3).getIntValue

      val teCalc: TransferEntropyCalculatorKraskov = new TransferEntropyCalculatorKraskov
      teCalc.setProperty("k", k.toString)
      teCalc.setProperty("delay", tau.toString)
      teCalc.initialise(1)
      teCalc.setObservations(X, Y)
      val res = teCalc.computeAverageLocalOfObservations()
      res.toLogoObject
    }
  }

  override def load(manager: PrimitiveManager) = {
    manager.addPrimitive("transfer-entropy", new TransferEntropy)
  }



}
